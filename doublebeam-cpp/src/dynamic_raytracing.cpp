#include "xtensor/xbroadcast.hpp"
#include <xtensor/xadapt.hpp>
#include <xtensor/xtensor.hpp>

#include "dynamic_raytracing.hpp"
#include "utils.hpp"


DynamicRayTracer::DynamicRayTracer(const VelocityModel& model) : kinematic(model), model(model) {}

using complex = std::complex<double>;

Beam DynamicRayTracer::trace_beam(state_type initial_state, double beam_width,
                                  double beam_frequency, const std::string& ray_code,
                                  double step_size, double max_step) {
    if (not model.in_model(initial_state[Index::Z])) {
        std::stringstream s;
        auto [top, bottom] = model.get_top_bottom();
        s << "Starting point of beam (" << initial_state[Index::X] << " " << initial_state[Index::Y]
          << " " << initial_state[Index::Z] << ") outside of model ( " << top << " " << bottom
          << ").";
        throw std::domain_error(s.str());
    }
    Beam beam(beam_width, beam_frequency);
    auto v0 = model.eval_at(initial_state[Index::Z]);
    xt::xtensor<complex, 3> P0{{{1j / v0, 0}, {0, 1j / v0}}};
    xt::xtensor<complex, 3> Q0{{{beam_frequency * beam_width * beam_width / v0, 0},
                                {0, beam_frequency * beam_width * beam_width / v0}}};
    auto ray = kinematic.trace_ray(initial_state, ray_code, step_size, max_step);
    for (auto& segment : ray) {
        // evaluate model at all points of the ray
        // TODO this was already done during ray tracing itself, maybe we can reuse the results
        std::vector<double> v;
        v.reserve(segment.data.size());
        for (auto& state : segment.data) {
            auto vel = model.eval_at(state[Index::Z]);
            v.push_back(vel * vel);
        }
        std::vector<double> travel_times(segment.data.size());
        std::transform(segment.data.begin(), segment.data.end(), travel_times.begin(),
                       [](const auto& el) { return el[Index::T]; });
        auto sigma_ = math::cumtrapz(v.begin(), v.end(), travel_times.begin(),
                                     travel_times.end(), 0.);
        auto sigma = xt::adapt(sigma_, {sigma_.size(), 1UL, 1UL});
        xt::xtensor<complex, 3> Q = Q0 + sigma * P0;
        xt::xtensor<complex, 3> P = xt::broadcast(P0, {sigma.size(), 2UL, 2UL});
        beam.segments.emplace_back(segment, P, Q);
    }

    return beam;
}
