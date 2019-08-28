#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <vector>
#include <array>
#include <type_traits>
#include <chrono>
#include "model.cpp"
#include "src/timing/timing.h"

typedef std::array<double, 7> state_type;


class KinematicRayTracingEquation{
    const VelocityModel& vm;

    double dvdz(double z){
        return 1 / vm[vm.layer_index(z)].gradient;
    }

public:
    explicit KinematicRayTracingEquation(const VelocityModel& model): vm(model) {}

    void operator()(const state_type& state, state_type& dxdt, const double d){
        auto [x, y, z, px, py, pz, T] = state;
        auto v = vm.eval_at(z);
        auto dxds = px * v;
        auto dyds = py * v;
        auto dzds = pz * v;
        auto dpzds = dvdz(z);
        auto dTds = 1. / v;
        dxdt[0] = dxds;
        dxdt[1] = dyds;
        dxdt[2] = dzds;
        dxdt[3] = 0.;
        dxdt[4] = 0.;
        dxdt[5] = dpzds;
        dxdt[6] = dTds;
    }
};


struct store_solver_state{
    std::vector<state_type>& states;
    std::vector<double>& times;

    store_solver_state(std::vector<state_type>& states, std::vector<double>& times): states(states), times(times){}

    void operator()(const state_type& x, double t){
        states.push_back(x);
        times.push_back(t);
    }
};


template <typename T, size_t N>
std::ostream& operator<<(std::ostream& os, std::array<T, N> a){
    os << "(";
    for (auto& i: a){
        os << i << " ";
    }
    os << ");";
    return os;
}


template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> v){
    for (auto& i: v){
        os << i << " ";
    }
    return os;
}


std::tuple<double, double, double> slowness_3D(double theta, double phi, double velocity) {
    auto px = 1. / velocity * sin(theta) * cos(phi);
    auto py = 1. / velocity * sin(theta) * sin(phi);
    auto pz = 1. / velocity * cos(theta);
    return {px, py, pz};
}

/*
template <typename T>
double radians(typename std::enable_if<std::is_integral<T>::value, T>::type degrees){
    constexpr T factor = M_PI / 180.;
    return degrees * factor;
}*/


template <typename T>
double radians(T degrees){
    constexpr double factor = M_PI / 180.;
    if constexpr (std::is_integral<T>::value){
        degrees = static_cast<double>(degrees);
    }
    return degrees * factor;
}


int main(){
    std::vector<Layer> layers{{0, 100, 1000, 1},
                              {100, 200, 1100, -1}};
    auto vm = VelocityModel(layers);
    const double startx = 0, starty = 0, startz = 0;
    const double velocity = vm.eval_at(startz);
    double theta = radians(50.);
    assert(theta != 0);
    double phi = radians(0.);
    const auto [px, py, pz] = slowness_3D(theta, phi, velocity);
    state_type initial_state{startx, starty, startz, px, py, pz, 0};
    std::vector<state_type> v;
    std::vector<double> t;

    auto trace = KinematicRayTracingEquation(vm);

    typedef boost::numeric::odeint::runge_kutta_cash_karp54< state_type > error_stepper_type;
    typedef boost::numeric::odeint::controlled_runge_kutta<error_stepper_type> controlled_stepper_type;
    controlled_stepper_type controlled_stepper;

    size_t steps = boost::numeric::odeint::integrate_adaptive(controlled_stepper, trace, initial_state, 0., 100., 1.,
            store_solver_state(v, t));
    std::cout << v << "\n";
    std::cout << steps << "\n";


    auto res = measure_runtime<8000>([&](){
        boost::numeric::odeint::integrate_adaptive(controlled_stepper, trace, initial_state, 0., 100., 1.,
                                                   store_solver_state(v, t));});
    std::cout << res << "\n";

    return 0;
}