#include "simulation.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace particle_life {
namespace {

constexpr double interaction_radius = 12.0;
constexpr double interaction_radius_sq = interaction_radius * interaction_radius;
constexpr double repulsion_radius = 2.0;
constexpr double force_scale = 0.018;
constexpr double damping = 0.91;
constexpr double max_speed = 1.35;
constexpr double time_step = 0.42;
constexpr double width = Simulation::width;
constexpr double height = Simulation::height;

}  // namespace

Simulation::Simulation(std::size_t particle_count, std::uint32_t seed)
    : particles_(particle_count), next_vx_(particle_count), next_vy_(particle_count) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> x_dist(0.0, width);
    std::uniform_real_distribution<double> y_dist(0.0, height);
    std::uniform_real_distribution<double> velocity_dist(-0.08, 0.08);
    std::uniform_real_distribution<double> affinity_dist(-1.0, 1.0);

    for (std::size_t source = 0; source < species_count; ++source) {
        for (std::size_t target = 0; target < species_count; ++target) {
            affinity_[source][target] = affinity_dist(rng);
        }
        affinity_[source][source] = 0.45 + 0.08 * static_cast<double>(source);
    }

    for (std::size_t i = 0; i < particle_count; ++i) {
        particles_[i] = {
            x_dist(rng),
            y_dist(rng),
            velocity_dist(rng),
            velocity_dist(rng),
            static_cast<std::uint8_t>(i % species_count),
        };
    }
}

double Simulation::wrapped_delta(double delta, double span) {
    if (delta > span * 0.5) {
        return delta - span;
    }
    if (delta < -span * 0.5) {
        return delta + span;
    }
    return delta;
}

double Simulation::wrap(double value, double span) {
    if (value >= span) {
        value -= span;
    } else if (value < 0.0) {
        value += span;
    }
    return value;
}

void Simulation::step() {
    // Uniform grid over the wrapped world. Cell sizes are at least interaction_radius, so every
    // neighbour within the radius lies in the 3x3 block of cells around a particle.
    constexpr std::size_t cols = static_cast<std::size_t>(width / interaction_radius);
    constexpr std::size_t rows = static_cast<std::size_t>(height / interaction_radius);
    constexpr double cell_w = width / static_cast<double>(cols);
    constexpr double cell_h = height / static_cast<double>(rows);
    static_assert(cols >= 3 && rows >= 3, "3x3 neighbourhood must not revisit cells");

    const std::size_t count = particles_.size();
    const std::size_t cell_count = cols * rows;

    // Counting-sort particles into cells (flat CSR layout).
    std::vector<std::size_t> particle_cell(count);
    std::vector<std::size_t> cell_start(cell_count + 1, 0);
    for (std::size_t i = 0; i < count; ++i) {
        const Particle& p = particles_[i];
        std::size_t cx = std::min(static_cast<std::size_t>(p.x / cell_w), cols - 1);
        std::size_t cy = std::min(static_cast<std::size_t>(p.y / cell_h), rows - 1);
        const std::size_t cell = cy * cols + cx;
        particle_cell[i] = cell;
        ++cell_start[cell + 1];
    }
    for (std::size_t c = 0; c < cell_count; ++c) {
        cell_start[c + 1] += cell_start[c];
    }
    std::vector<std::size_t> cell_items(count);
    {
        std::vector<std::size_t> fill(cell_start.begin(), cell_start.end() - 1);
        for (std::size_t i = 0; i < count; ++i) {
            cell_items[fill[particle_cell[i]]++] = i;
        }
    }

    std::vector<std::size_t> candidates;
    for (std::size_t i = 0; i < count; ++i) {
        const Particle& particle = particles_[i];
        double force_x = 0.0;
        double force_y = 0.0;

        // Gather candidates from the 3x3 wrapped neighbourhood, then sort so that the
        // accumulation order matches the original ascending-j loop.
        candidates.clear();
        const std::size_t cell = particle_cell[i];
        const std::size_t cx = cell % cols;
        const std::size_t cy = cell / cols;
        for (std::size_t oy = 0; oy < 3; ++oy) {
            const std::size_t ny = (cy + rows + oy - 1) % rows;
            for (std::size_t ox = 0; ox < 3; ++ox) {
                const std::size_t nx = (cx + cols + ox - 1) % cols;
                const std::size_t nc = ny * cols + nx;
                for (std::size_t k = cell_start[nc]; k < cell_start[nc + 1]; ++k) {
                    if (cell_items[k] != i) {
                        candidates.push_back(cell_items[k]);
                    }
                }
            }
        }
        std::sort(candidates.begin(), candidates.end());

        for (const std::size_t j : candidates) {
            const Particle& other = particles_[j];
            const double dx = wrapped_delta(other.x - particle.x, width);
            const double dy = wrapped_delta(other.y - particle.y, height);
            const double distance_sq = dx * dx + dy * dy;

            if (distance_sq <= 1.0e-12 || distance_sq >= interaction_radius_sq) {
                continue;
            }

            const double distance = std::sqrt(distance_sq);
            double strength;
            if (distance < repulsion_radius) {
                strength = distance / repulsion_radius - 1.0;
            } else {
                const double envelope = 1.0 - std::abs(2.0 * distance / interaction_radius - 1.0);
                strength = affinity_[particle.species][other.species] * envelope;
            }

            const double scaled = force_scale * strength / distance;
            force_x += dx * scaled;
            force_y += dy * scaled;
        }

        double vx = (particle.vx + force_x) * damping;
        double vy = (particle.vy + force_y) * damping;
        const double speed_sq = vx * vx + vy * vy;
        if (speed_sq > max_speed * max_speed) {
            const double scale = max_speed / std::sqrt(speed_sq);
            vx *= scale;
            vy *= scale;
        }

        next_vx_[i] = vx;
        next_vy_[i] = vy;
    }

    for (std::size_t i = 0; i < particles_.size(); ++i) {
        Particle& particle = particles_[i];
        particle.vx = next_vx_[i];
        particle.vy = next_vy_[i];
        particle.x = wrap(particle.x + particle.vx * time_step, width);
        particle.y = wrap(particle.y + particle.vy * time_step, height);
    }
}

double Simulation::checksum() const {
    long double result = 0.0;
    for (std::size_t i = 0; i < particles_.size(); ++i) {
        const Particle& particle = particles_[i];
        const long double weight = static_cast<long double>((i % 97) + 1);
        result += weight * (particle.x * 0.7 + particle.y * 1.3 + particle.vx * 2.1 +
                            particle.vy * 2.9 + static_cast<double>(particle.species) * 0.17);
    }
    return static_cast<double>(result);
}

}  // namespace particle_life
