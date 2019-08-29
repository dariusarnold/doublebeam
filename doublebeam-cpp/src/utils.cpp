#include "utils.h"



std::tuple<double, double, double> seismo::slowness_3D(double theta, double phi, double velocity){
    auto px = 1. / velocity * sin(theta) * cos(phi);
    auto py = 1. / velocity * sin(theta) * sin(phi);
    auto pz = 1. / velocity * cos(theta);
    return {px, py, pz};
}