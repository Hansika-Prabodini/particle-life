#include "simulation.hpp"

#include <cmath>
#include <iostream>

namespace {

bool close(double left, double right, double tolerance = 1.0e-12) {
    return std::abs(left - right) <= tolerance;
}

}  // namespace

int main() {
    particle_life::Simulation single_particle(1, 7);
    auto& particle = const_cast<particle_life::Particle&>(single_particle.particles().front());
    particle.x = 3.0 * particle_life::Simulation::width - 0.1;
    particle.y = -2.0 * particle_life::Simulation::height + 0.1;
    particle.vx = 0.0;
    particle.vy = 0.0;
    single_particle.step();

    constexpr double expected_x = particle_life::Simulation::width - 0.1;
    constexpr double expected_y = 0.1;
    if (!close(particle.x, expected_x) || !close(particle.y, expected_y)) {
        std::cerr << "positions crossing multiple world spans were not wrapped\n";
        return 1;
    }

    particle_life::Simulation first(128, 42);
    particle_life::Simulation second(128, 42);

    for (int step = 0; step < 20; ++step) {
        first.step();
        second.step();
    }

    const auto& first_particles = first.particles();
    const auto& second_particles = second.particles();
    if (first_particles.size() != 128 || second_particles.size() != 128) {
        std::cerr << "particle count changed\n";
        return 1;
    }

    for (std::size_t i = 0; i < first_particles.size(); ++i) {
        const auto& left = first_particles[i];
        const auto& right = second_particles[i];
        if (!close(left.x, right.x) || !close(left.y, right.y) || !close(left.vx, right.vx) ||
            !close(left.vy, right.vy) || left.species != right.species) {
            std::cerr << "simulation is not deterministic at particle " << i << '\n';
            return 1;
        }
        if (left.x < 0.0 || left.x >= particle_life::Simulation::width || left.y < 0.0 ||
            left.y >= particle_life::Simulation::height || !std::isfinite(left.vx) ||
            !std::isfinite(left.vy)) {
            std::cerr << "particle escaped simulation bounds at index " << i << '\n';
            return 1;
        }
    }

    particle_life::Simulation wrapping(1, 7);
    for (int step = 0; step < 200; ++step) {
        wrapping.step();
    }
    const auto& wrapped_particle = wrapping.particles().front();
    if (wrapped_particle.x < 0.0 || wrapped_particle.x >= particle_life::Simulation::width ||
        wrapped_particle.y < 0.0 || wrapped_particle.y >= particle_life::Simulation::height) {
        std::cerr << "particle did not wrap after a large displacement\n";
        return 1;
    }

    const double result = first.checksum();
    constexpr double expected_checksum = 610659.60219911579;
    if (!std::isfinite(result) || !close(result, expected_checksum, 1.0e-6)) {
        std::cerr.precision(17);
        std::cerr << "semantic checksum changed: expected " << expected_checksum << ", got " << result
                  << '\n';
        return 1;
    }

    std::cout.precision(17);
    std::cout << "deterministic checksum: " << result << '\n';
    return 0;
}
