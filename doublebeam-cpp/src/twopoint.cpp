#include "twopoint.hpp"


TwoPointRayTracing::TwoPointRayTracing(VelocityModel& velocity_model) : model(velocity_model) {}

std::string stringify(position_t pos) {
    auto [x, y, z] = pos;
    return "(" + std::to_string(x) + " " + std::to_string(y) + " " + std::to_string(z) + ")";
}

slowness_t TwoPointRayTracing::trace(position_t source, position_t receiver, double accuracy) {
    auto [top, bottom] = model.get_top_bottom();
    auto [source_x, source_y, source_z] = source;
    auto [receiver_x, receiver_y, receiver_z] = receiver;
    if (not model.in_model(source_z)) {
        throw std::domain_error("Source outside of model " + stringify(source));
    }
    if (not model.in_model(receiver_z)) {
        throw std::domain_error("Receiver outside of model " + stringify(receiver));
    }


    return slowness_t();
}
