#include <fmt/format.h>
#include <fmt/ostream.h>

#include <raytracing_types.hpp>


std::ostream& operator<<(std::ostream& os, Position position) {
    return os << fmt::format("Position({} m, {} m, {} m)", position.x.get(), position.y.get(),
                             position.z.get());
}

std::ostream& operator<<(std::ostream& os, Slowness slowness) {
    return os << fmt::format("Slowness({} s/m, {} s/m, {} s/m)", slowness.px.get(),
                             slowness.py.get(), slowness.pz.get());
}

std::ostream& operator<<(std::ostream& os, TravelTime travel_time) {
    return os << fmt::format("Traveltime({} s)", travel_time.time.get());
}

std::ostream& operator<<(std::ostream& os, Arclength arclength) {
    return os << fmt::format("Arclength({} m)", arclength.length.get());
}

std::ostream& operator<<(std::ostream& os, RayState ray_state) {
    return os << fmt::format("RayState({}, {}, {}, {})", ray_state.position, ray_state.slowness,
                             ray_state.travel_time, ray_state.arclength);
}
