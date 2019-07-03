import enum
from collections import namedtuple
from math import sin, cos, asin, acos, radians, degrees, copysign, sqrt
from typing import Tuple, Callable, List, Union

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate
import scipy.misc
from scipy.integrate import solve_ivp
from scipy.misc import derivative

from doublebeam.core.models import VelocityModel3D, LinearVelocityLayer


def cartesian_to_ray_s(x, z, xm, _theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    s, see eq. 14a in Hill1990 Gaussian beam migration"""
    return (x - xm) * sin(_theta) - z * cos(_theta)


def cartesian_to_ray_n(x, z, xm, theta):
    """Conversion of 2D cartesian coordinates x, z to ray coordinate
    n, see eq. 14b in Hill1990 Gaussian beam migration"""
    return (x - xm) * cos(theta) + z * sin(theta)


def dvx():
    """Derivative of velocity after x. 0 For 1D model"""
    return 0.


def dvz(velocity_model, z, delta=0.0001):
    """Derivative of velocity after z"""
    if z < delta:
        # special case: evaluating derivative would eval model at negative depth
        return 0.
    return scipy.misc.derivative(velocity_model.eval_at, z, delta)


def horizontal_slowness(v, theta):
    return sin(theta) / v


def vertical_slowness(v, theta):
    return cos(theta) / v


class Ray2D:

    def __init__(self, start_x: float, start_z: float, theta: float):
        """
        :param start_x: x coordinate of start point of ray in m
        :param start_z: z coordinate of start point of ray in m
        :param theta: angle of ray against downgoing vertical (z axis) at the
        surface in rad
        """
        self.x0 = start_x
        self.z0 = start_z
        self.theta = theta
        # Stores x, y coordinates of ray path
        self.path = None  # type: Tuple[np.array, np.array]
        self.layer_boundaries_crossed_depths = [-99]

    def get_last_boundary_crossed_depth(self):
        try:
            return self.layer_boundaries_crossed_depths[-1]
        except IndexError:
            return None


class Ray3D:

    def __init__(self, start_x: float, start_y: float, start_z: float, theta: float, phi: float):
        """
        :param start_x: x coordinate of start point of ray in m
        :param start_y: y coordinate of start point of ray in m
        :param start_z: z coordinate of start point of ray in m
        :param theta: Angle against downgoing vertical axis (z) at start
        point in rad, increasing upwards. 0 <= theta <= pi
        :param phi: Angle against x axis at start point in rad, with increasing
        angle towards the y axis
        0 <= phi <= 2*pi
        """
        self.start = np.array((start_x, start_y, start_z))
        self.theta = theta
        self.phi = phi
        # These will be set after the ray is traced
        self.path: List[np.ndarray] = []
        self.slowness: List[np.ndarray] = []
        self.travel_time: List[np.ndarray] = []

    @property
    def last_point(self) -> np.ndarray:
        """
        Return last point of the rays path. If the ray hasn't been traced yet,
        return the starting point of the ray.
        """
        try:
            return self.path[-1][-1]
        except IndexError:
            return self.start

    @property
    def last_slowness(self) -> Union[np.ndarray, None]:
        """
        Return last slowness value. If the ray hasn't been traced yet, return
        None.
        """
        try:
            return self.slowness[-1][-1]
        except IndexError:
            return None

    @property
    def last_time(self) -> float:
        """
        Return last travel time of the ray.
        :return:
        """
        try:
            return self.travel_time[-1][-1]
        except IndexError:
            return 0.

    @property
    def direction(self) -> str:
        """
        Return direction of the ray along the z-axis. A downgoing ray will
        return "down", while an upgoing ray will return "up". The special case
        of a horizontal ray (vertical slowness pz = 0) will return "horizontal".
        """
        try:
            pz = self.last_slowness[Index.Z]
        except TypeError:
            # if None is returned (ray hasn't been traced yet), use starting
            # angle against vertical axis
            pz = np.pi/2 - self.theta
        if pz < 0:
            return "up"
        elif pz > 0:
            return "down"
        return "horizontal"


def calc_initial_slowness3D(ray: Ray3D, v0: float) -> Tuple[float, float, float]:
    """
    Calculate initial vertical and horizontal slowness for a ray.
    For geometric definitions see chapter 3.2.1 in Cerveny - Seismic ray theory
    :param ray: Ray instance
    :param v0: Velocity at ray starting point
    :return: Tuple of slowness values
    """
    px = 1/v0 * sin(ray.theta) * cos(ray.phi)
    py = 1/v0 * sin(ray.theta) * sin(ray.phi)
    pz = 1/v0 * cos(ray.theta)
    return px, py, pz


def length(vector: np.ndarray) -> float:
    """Calculate length of vector"""
    return (vector @ vector)**0.5


def angle(vector1: np.ndarray, vector2: np.ndarray) -> float:
    """Calculate angle between two vectors in rad"""
    return acos((vector1 @ vector2) / (length(vector1) * length(vector2)))


def critical_angle(v1: float, v2: float) -> float:
    """
    Use Snells law to calculate the critical angle at an interface
    :param v1: velocity before interface in m/s
    :param v2: velocity after interface in m/s
    :return: critical angle in rad
    """
    if v1 < v2:
        return asin(v1/v2)
    return np.pi


def snells_law_(px: float, pz: float, z: float, z_prev: float, velocity_model: VelocityModel3D, wave_type: str = "T"):
    """
    Compute new slowness vector using Snells law
    :param px: x component of slowness vector before interface
    :param pz: z component of slowness vector before interface
    :param z: depth after interface
    :param z_prev: depth before interface
    :param wave_type: Select for which wave type (reflected/transmitted) the
    slowness should be computed. "T" for transmitted, "R" for reflected.
    :return: slowness vector (x, z) after interface
    """
    # create slowness vector p
    # TODO keep p as a np.ndarray everywhere
    p = np.array((px, pz))
    # normal vector in 1D model will always be vertical
    # n should be oriented to the side where the transmitted wave propagates
    # else the minus/plus relation for transmitted/reflected waves isn't valid
    n = np.array((0, copysign(1, pz)))
    # look at angle to determine critical incidence
    angle_in = angle(p, n)
    print("Angle in", degrees(angle_in))
    eps = copysign(1, np.dot(p, n))
    v = velocity_model.eval_at(z)
    v_prev = velocity_model.eval_at(z_prev)
    angle_crit = critical_angle(v_prev, v)
    print("Critical angle", degrees(angle_crit))
    if angle_in > angle_crit:
        print("Total reflection")
        p[1] *= -1
        return p
    #angle_out = asin(sin(angle_in) * v/v_prev)
    #print("Angle out", degrees(angle_out))
    pn = np.dot(p, n)
    minusplus = -1 if wave_type == "T" else 1
    # TODO not stable when using reflected wave
    # calculate slowness vector after interface using eq. 2.4.70 from
    # Cerveny - Seismic ray theory
    p_new = p - n * (pn + minusplus * eps * sqrt(v**-2 - v_prev**-2 + pn**2))
    return p_new


def snells_law(p: np.ndarray, v_before: float, v_after: float, wave_type: str = "T") -> np.ndarray:
    """
    Calculate initial slowness of a wave transmitted across an interface.
    Eq. (2.4.70) in Cerveny - Seismic ray theory (2001).
    :param p: Slowness vector before interface
    :param v_before: Velocity before the interface
    :param v_after: Velocity after the interface
    :param wave_type: Select for which wave type (reflected/transmitted) the
    slowness should be computed. "T" for transmitted, "R" for reflected.
    """
    # n is normal vector of the interface. n should be oriented to the side
    # the transmitted wave propagates, else the minus/plus relation for
    # transmitted/reflected waves isn't valid.
    n = np.array((0, 0, copysign(1, p[2])))
    pn = p.dot(n)
    eps = copysign(1, pn)
    minusplus = -1 if wave_type == "T" else 1
    return p - (pn + minusplus * eps * (v_after**-2 - v_before**-2 + pn**2)**0.5) * n


def plot_ray(x1: np.ndarray, x2: np.ndarray):
    """Plot two coordinates of a aray"""
    plt.plot(x1, x2, label="Ray path")
    ax = plt.gca()
    # invert y axis so positive depth values are shown downwards
    ax.invert_yaxis()
    # set aspect ratio to equal so angles stay true
    ax.set_aspect("equal")
    plt.xlabel("x1 (m)")
    plt.ylabel("x2 (m)")
    plt.legend()
    plt.show()


def plot_ray_in_model(ray: Ray3D, velocity_model: VelocityModel3D):
    fig, ax = plt.subplots()
    ax.invert_yaxis()
    ax.set_ylabel("Depth (m)")
    ax.set_xlabel("Offset along x (m)")
    for depth in velocity_model.interface_depths:
        ax.axhline(depth, color="k")
    for segment in ray.path:
        x, z = segment.T[Index.X], segment.T[Index.Z]
        ax.plot(x, z)
    plt.show()


class IVPResultStatus(enum.IntEnum):
    """Enum class to make the int status returned from scipys solve_ivp more readable"""
    FAILED = -1
    END_REACHED = 0
    TERMINATION_EVENT = 1


ODEState = namedtuple("ODEState", ["x", "z", "px", "pz"])

IVPEventFunction = Callable[[float, ODEState], float]

class RayTracer2D:
    """
    Class for ray tracing in a 2D velocity model.
    """

    def __init__(self, velocity_model: VelocityModel3D):
        """
        :param velocity_model: Velocity model to use for ray tracing
        """
        self._velocity_model = velocity_model
        # Generate a function which has a zero crossing at the boundary depth
        # for all boundary depths in the models to apply Snells law
        # skip first interface depth since its the surface at z = 0 and skip
        # last interface since model stops there
        #  TODO add condition that stops integration once model bottom is reached
        crossings : List[IVPEventFunction] = [lambda s, y: y[1] - depth
                                              for depth in velocity_model.interface_depths[1:-1]]
        # set to True to stop integration at the boundary so we can apply Snells law
        for f in crossings:
            f.terminal = True
            f.direction = 1
        # return z coordinate which has a natural zero crossing at the surface
        surfaced: IVPEventFunction = lambda s, y: y[1]
        # set True to stop integration once the ray reaches the surface
        surfaced.terminal = True
        self._events = [*crossings, surfaced]

    def _trace(self, s: float, y: ODEState) -> ODEState:
        """
        Standard raypath equations in 2D from Hill1990 Gaussian beam migration eq. 2a-2d
        :param s: Current value of integration variable s (arclength along ray)
        :param y: List of values for x, y coordinate, horizontal and vertical slowness
        :return:
        """
        x, z, px, pz = y
        try:
            v = self._velocity_model.eval_at(z)
        except LookupError:
            # evaluating at negative depth means surface is reached
            # TODO more clear/expressive way to stop integration
            return ODEState(0, 0, 0, 0)
        dxds = v * px
        dzds = v * pz
        # TODO simplify dvx by replacing it with zero for 2D case
        dpxds = -1 * v**-2 * dvx()
        dpzds = -1 * v**-2 * dvz(self._velocity_model, z)
        dydt = ODEState(dxds, dzds, dpxds, dpzds)
        return dydt

    def ray_trace(self, ray: Ray2D, s_end: float, ds: float=0.01) -> Ray2D:
        v0 = self._velocity_model.eval_at(ray.z0)
        px0 = horizontal_slowness(v0, ray.theta)
        pz0 = vertical_slowness(v0, ray.theta)
        min_float_step = np.finfo(float).eps
        initial_state = ODEState(ray.x0, ray.z0+min_float_step, px0, pz0)
        result: scipy.integrate._ivp.ivp.OdeResult = solve_ivp(self._trace, [0, s_end], initial_state, max_step=ds,
                           events=self._events)
        x_values, z_values = [], []
        while result.status == IVPResultStatus.TERMINATION_EVENT:
            # events are returned in the order as passed to solve_ivp
            crossing_events, surface_events = result.t_events
            # stop integration when surface was reached
            if surface_events.size > 0:
                break
            s_event = crossing_events[0]
            x_, z_, px_, pz_ = result.y
            x_values.append(x_)
            z_values.append(z_)
            px, pz = snells_law(px_[-1], pz_[-1], z_[-1], z_[-2], self._velocity_model)
            min_float_step = copysign(min_float_step, pz_[-1])
            # move z behind the interface just passed by the ray so the event wont
            # trigger again. Increase z (positive step) when ray goes down, decrease
            # z (negative step) when ray goes up
            #z_[-1] += copysign(z_[-1], 50000)
            result = solve_ivp(self._trace, [s_event, s_end],
                               ODEState(x_[-1], z_[-1], px, pz),
                               max_step=ds, events=self._events)
        x_values.append(result.y[0])
        z_values.append(result.y[1])
        x_values = np.concatenate(x_values)
        z_values = np.concatenate(z_values)
        ray.path = (x_values, z_values)
        return ray


ODEState3D = namedtuple("ODEState3D", ["x", "y", "z", "px", "py", "pz", "T"])


class Index(enum.IntEnum):
    X = 0
    Y = 1
    Z = 2


class NumericRayTracer3D:

    def __init__(self, velocity_model: VelocityModel3D):
        self.velocity_model = velocity_model
        self.layer = None

    @staticmethod
    def _velocity(layer: LinearVelocityLayer, x: float, y: float, z: float) -> float:
        """
        Evaluate velocity at a given depth in a layer with linearly varying velocity
        :param layer: The layer to evaluate
        :param depth: Depth in m at which to evaluate the property
        """
        return layer["intercept"] + layer["gradient"] * z

    def _trace(self, s: float, y: ODEState3D) -> ODEState3D:
        """
        Implement ray tracing system (3.1.10) from Cerveny - Seismic ray theory
        (2001).
        :param s: Arclength along ray
        :param y: Previous state
        :return: New state
        """
        x, y, z, px, py, pz, T = y
        v = self._velocity(self.layer, x, y, z)
        dxds = px * v
        dyds = py * v
        dzds = pz * v
        dpxds = derivative((lambda x: 1. / self._velocity(self.layer, x, y, z)), x, dx=0.0001)
        dpyds = derivative((lambda y: 1. / self._velocity(self.layer, x, y, z)), y, dx=0.0001)
        dpzds = derivative((lambda z: 1. / self._velocity(self.layer, x, y, z)), z, dx=0.0001)
        dTds = 1. / v
        return ODEState3D(dxds, dyds, dzds, dpxds, dpyds, dpzds, dTds)

    def trace_layer(self, ray: Ray3D, initial_slownesses: Tuple[float, float, float]) -> Ray3D:
        """
        Trace a ray through a single layer of the model
        :param ray: Ray to trace
        :param initial_slownesses: Initial value of px, py, pz (slowness along
        the corresponding axis) of the ray in s/m.
        :return: Ray with its parameters set
        """
        out_of_layer_events = [lambda s, y: y[2] - depth for depth in
                               (self.layer["top_depth"], self.layer["bot_depth"])]
        for func in out_of_layer_events:
            func.terminal = True
        # TODO change signature of _velocity to accept np array instead of three floats
        initial_state = ODEState3D(*ray.last_point, *initial_slownesses, ray.last_time)
        result: scipy.integrate._ivp.ivp.OdeResult = solve_ivp(self._trace, (0, np.inf), initial_state, max_step=1, events=out_of_layer_events)
        x, y, z, px, py, pz, t = result.y
        ray.path.append(np.vstack((x, y, z)).T)
        ray.slowness.append(np.vstack((px, py, pz)).T)
        ray.travel_time.append(t)
        return ray

    def trace_stack(self, ray: Ray3D, ray_code: str = None) -> Ray3D:
        """
        Trace ray through a stack of layers. The ray type at an interface is
        chosen by the ray code.
        :param ray: Ray to trace through the model
        :param ray_code: Specifies which ray (Transmitted/Reflected) to follow
        at an interface in the model. "T" stands for the transmitted ray, "R"
        for the reflected ray. If not given, the ray will only be traced through
        the layer in which its starting point resides.
        :return: Ray with parameters set
        """
        index = self.velocity_model.layer_index(ray.start[Index.Z])
        self.layer = self.velocity_model[index]
        initial_slownesses = calc_initial_slowness3D(ray, self._velocity(self.layer, *ray.last_point))
        ray = self.trace_layer(ray, initial_slownesses)
        if ray_code is None:
            return ray
        for wave_type in ray_code:
            v_top, v_bottom = self.velocity_model.interface_velocities(ray.last_point[Index.Z])
            index = self.velocity_model.layer_index(ray.last_point[Index.Z])
            self.layer = self.velocity_model[index]
            if ray.direction == "down":
                new_p = snells_law(ray.last_slowness, v_top, v_bottom, wave_type)
            else:
                new_p = snells_law(ray.last_slowness, v_bottom, v_top, wave_type)
            ray = self.trace_layer(ray, new_p)
        return ray


def main():
    layers = [(0, 100, 1800, 4), (100, 200, 2400, 0), (200, 300, 2400, 1),
              (300, 400, 2700, 0), (400, 500, 2250, 1.5)]
    vm = VelocityModel3D(layers)
    ray = Ray3D(0, 0, 0, radians(20), radians(0))
    nrt = NumericRayTracer3D(vm)
    ray = nrt.trace_stack(ray, "TTTTRTTTT")
    plot_ray_in_model(ray, vm)


if __name__ == '__main__':
    main()
