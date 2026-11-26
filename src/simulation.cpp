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
    const std::size_t n = particles_.size();

    // Uniform toroidal grid with cell size >= interaction_radius.
    int cols = static_cast<int>(std::floor(width / interaction_radius));
    int rows = static_cast<int>(std::floor(height / interaction_radius));
    if (cols < 1) cols = 1;
    if (rows < 1) rows = 1;
    const double cell_w = width / static_cast<double>(cols);
    const double cell_h = height / static_cast<double>(rows);
    const std::size_t cell_count = static_cast<std::size_t>(cols) * static_cast<std::size_t>(rows);

    cell_of_.resize(n);
    cell_start_.assign(cell_count + 1, 0);
    cell_items_.resize(n);
    candidates_.reserve(n);

    for (std::size_t i = 0; i < n; ++i) {
        int cx = static_cast<int>(particles_[i].x / cell_w);
        int cy = static_cast<int>(particles_[i].y / cell_h);
        cx = std::min(std::max(cx, 0), cols - 1);
        cy = std::min(std::max(cy, 0), rows - 1);
        const std::uint32_t c = static_cast<std::uint32_t>(cy * cols + cx);
        cell_of_[i] = c;
        ++cell_start_[c + 1];
    }
    for (std::size_t c = 0; c < cell_count; ++c) {
        cell_start_[c + 1] += cell_start_[c];
    }
    {
        // Counting sort: indices within each cell are stored in ascending order.
        std::vector<std::uint32_t>& fill = candidates_;
        fill.assign(cell_start_.begin(), cell_start_.end() - 1);
        for (std::size_t i = 0; i < n; ++i) {
            cell_items_[fill[cell_of_[i]]++] = static_cast<std::uint32_t>(i);
        }
    }

    // Deduplicated neighbour offsets per axis.
    int col_offsets[3];
    int row_offsets[3];
    int ncol_off = cols >= 3 ? 3 : cols;
    int nrow_off = rows >= 3 ? 3 : rows;

    for (std::size_t i = 0; i < n; ++i) {
        const Particle& particle = particles_[i];
        double force_x = 0.0;
        double force_y = 0.0;

        const int ci = static_cast<int>(cell_of_[i]);
        const int cx = ci % cols;
        const int cy = ci / cols;
        if (cols >= 3) {
            col_offsets[0] = (cx + cols - 1) % cols;
            col_offsets[1] = cx;
            col_offsets[2] = (cx + 1) % cols;
        } else {
            for (int k = 0; k < cols; ++k) col_offsets[k] = k;
        }
        if (rows >= 3) {
            row_offsets[0] = (cy + rows - 1) % rows;
            row_offsets[1] = cy;
            row_offsets[2] = (cy + 1) % rows;
        } else {
            for (int k = 0; k < rows; ++k) row_offsets[k] = k;
        }

        candidates_.clear();
        for (int r = 0; r < nrow_off; ++r) {
            for (int c = 0; c < ncol_off; ++c) {
                const std::size_t cell =
                    static_cast<std::size_t>(row_offsets[r]) * static_cast<std::size_t>(cols) +
                    static_cast<std::size_t>(col_offsets[c]);
                candidates_.insert(candidates_.end(), cell_items_.begin() + cell_start_[cell],
                                   cell_items_.begin() + cell_start_[cell + 1]);
            }
        }
        std::sort(candidates_.begin(), candidates_.end());

        for (const std::uint32_t j32 : candidates_) {
            const std::size_t j = j32;
            if (i == j) {
                continue;
            }

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
