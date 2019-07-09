"""
Two point ray tracing as described by
A fast and robust two-point ray tracing method in layered media with constant or
linearly varying layer velocity (Xinding Fang and Xiaofei Chen, 2019)
"""
import time
from math import sqrt

import numpy as np

from doublebeam.core.raytracing.raytracing import VelocityModel3D
from doublebeam.core.common import Index


class TwoPointRayTracing:

    def __init__(self, velocity_model: VelocityModel3D):
        self._velocity_model: VelocityModel3D = velocity_model
        self.a = velocity_model.gradients
        self.b = velocity_model.intercepts
        self.z = velocity_model.interface_depths
        self.n = len(velocity_model.gradients)

    def _epsilon(self, k: int, s: int) -> float:
        """
        Eq. A5
        :param k: Index of current layer
        :param s: Index of source layer
        """
        if k == 0:
            s += 1
            return (self.a[s] * self.z[s-1] + self.b[k])**2
        else:
            return (self.a[k] * self.z[k-1] + self.b[k])**2

    def _omega(self, k: int, s: int) -> float:
        """
        Eq. A6
        :param k: Index of current layer
        :param s: Index of source layer
        """
        if k == 0:
            s += 1
            return (self.a[s] * self.zs + self.b[s])**2
        else:
            return (self.a[k] * self.z[k] + self.b[k])**2

    def _h(self, k: int, s: int) -> float:
        """
        Eq. A7
        :param k: Index of current layer
        :param s: Index of source layer
        """
        if k == 0:
            s += 1
            return (self.a[s] * self.z[s-1] + self.b[s]) * (self.zs - self.z[s-1])
        else:
            return (self.a[k] * self.z[k-1] + self.b[k]) * (self.z[k] - self.z[k-1])

    def _mu(self, k: int, s: int, wave_type: str = "direct") -> float:
        """
        Eq. A8
        :param k: Index of current layer
        :param s: Index of layer that contains the source
        :param wave_type: "direct" or "reflected"
        """
        if k == 0:
            s += 1
            return 1 - self._mu(k=s, s=s, wave_type=wave_type)
        elif 1 <= k <= s-1:
            return 1.
        elif s <= k <= self.n:
            return 0. if wave_type == "direct" else 2.

    def _delta_a(self, k: int) -> float:
        """
        Eq. A9
        :param k: Index of current layer
        """
        return 1 if self.a[k] != 0. else 0.

    def _mu_tilde(self, k: int, s: int) -> float:
        """
        Eq. A10
        :param k: Index of current layer
        :param s: Index of source layer
        """
        if self.a[k] == 0:
            raise ValueError("a_k is zero")
        return self._mu(k, s) * self.v_M / self.a[k]

    def _h_tilde(self, k: int, s: int) -> float:
        """
        Eq. A11
        :param k: Index of current layer
        :param s: Index of source layer
        """
        # TODO add wave type parameter to _mu call
        return self._mu(k, s) * self._h(k, s) / self.v_M

    def _epsilon_tilde(self, k: int, s: int) -> float:
        """
        Eq. A12
        :param k: Index of current layer
        :param s: Index of source layer
        """
        return 1 - self._epsilon(k, s) / self.v_M**2

    def _omega_tilde(self, k: int, s: int) -> float:
        """
        Eq. A13
        :param k: Index of current layer
        :param s: Index of source layer
        """
        return 1 - self._omega(k, s) / self.v_M**2

    def _delta_epsilon(self, k: int, s: int) -> float:
        """
        Eq. C18
        :param k: Index of current layer
        :param s: Index of source layer
        """
        epsilon_tilde = self._epsilon_tilde(k, s)
        return 1 if epsilon_tilde != 0 else 0

    def _delta_omega(self, k: int, s: int) -> float:
        """
        Eq. C19
        :param k: Index of current layer
        :param s: Index of source layer
        """
        omega_tilde = self._omega_tilde(k, s)
        return 1 if omega_tilde != 0 else 0

    def _q_A(self) -> float:
        """
        Eq. C11
        """
        return 1 / (self.v_A * sqrt(self.v_A**-2 - self.v_M**-2))

    def _d0(self, s: int) -> float:
        """
        Eq. C12
        :param s: Index of source layer
        """
        def d0_iter(k: int) -> float:
            """
            Inner part of sum from eq. C12
            :param k: Index of current layer
            """
            return ((self._mu_tilde(k, s)
                    * (sqrt(self._epsilon_tilde(k, s) - self._q_A()**-2)
                       - sqrt(self._omega_tilde(k, s) - self._q_A()**-2)) if self._delta_a(k) else 0.)
                    + (self._h_tilde(k, s) / sqrt(self._epsilon_tilde(k, s) - self._q_A()**-2) if (1 - self._delta_a(k)) else 0.))

        return (sum(d0_iter(k) for k in range(0, self.n))
                + self._mu_tilde(self.n, s) * sqrt(self._epsilon_tilde(self.n, s)
                                                   - self._q_A()**-2))

    def _d1(self, s: int) -> float:
        """
        Eq. C13
        :param s: Index of source layer
        """
        def d1_iter(k: int) -> float:
            """
            Inner part of sum from eq. C13
            :param k: Index of current layer
            """
            return ((0.5 * self._mu_tilde(k, s) * (self._epsilon_tilde(k, s) - self._omega_tilde(k, s)) if self._delta_a(k) else 0.)
                    + (1 - self._delta_a(k)) * self._h_tilde(k, s))

        return sum(d1_iter(k) for k in range(0, self.n+1))

    def _c0(self, s: int) -> float:
        """
        Eq. C14
        :param s: Index of source layer
        """
        def c0_iter(k: int) -> float:
            """
            Inner part of sum from eq. C14
            :param k: Index of current layer
            """
            return ((self._mu_tilde(k, s) * ((sqrt(self._epsilon_tilde(k, s)) if self._delta_epsilon(k, s) else 0.)
                                             - (sqrt(self._omega_tilde(k, s)) if self._delta_omega(k, s) else 0.)) if self._delta_a(k) else 0.)
                    + (self._h_tilde(k, s) / sqrt(self._epsilon_tilde(k, s)) if (1 - self._delta_a(k)) * self._delta_epsilon(k, s) else 0.))

        return sum(c0_iter(k) for k in range(0, self.n+1))

    def _c1(self, s: int) -> float:
        """
        Eq. C15
        :param s: Index of source layer
        """
        def c1_iter(k: int) -> float:
            """
            Inner part of sum from eq. C15
            :param k: Index of current layer
            """
            return self._h_tilde(k, s) if (1 - self._delta_a(k)) * (1 - self._delta_epsilon(k, s)) else 0.

        return sum(c1_iter(k) for k in range(0, self.n+1))

    def _c_minus1(self, s: int) -> float:
        """
        Eq. C16
        :param s: Index of source layer
        """
        def c_minus1_iter(k: int) -> float:
            """
            Inner part of sum from eq. C16
            :param k: Index of current layer
            """
            return self._mu_tilde(k, s) if self._delta_a(k) * (self._delta_omega(k, s) - self._delta_epsilon(k, s)) else 0.

        return sum(c_minus1_iter(k) for k in range(0, self.n+1))

    def _c_minus2(self, s: int) -> float:
        """
        Eq. C17
        :param s: Index of source layer
        """
        def c_minus2_iter(k: int) -> float:
            """
            Inner part of sum from eq. C17
            :param k: Index of current layer
            """
            return ((0.5 * self._mu_tilde(k, s)
                    * ((1 / sqrt(self._epsilon_tilde(k, s)) if self._delta_epsilon(k, s) else 0.)
                       - (1 / sqrt(self._omega_tilde(k, s)) if self._delta_omega(k, s) else 0.)) if self._delta_a(k) else 0.)
                    - (0.5 * self._h_tilde(k, s) / self._epsilon_tilde(k, s)**1.5 if (1 - self._delta_a(k)) * self._delta_epsilon(k, s) else 0.))

        return sum(c_minus2_iter(k) for k in range(0, self.n+1))

    def _alpha1(self, s: int) -> float:
        """
        Eq. C2
        :param s: Index of source layer
        """
        return self._d1(s)

    def _alpha2(self, s: int) -> float:
        """
        Eq. C3
        :param s: Index of source layer
        """
        c0 = self._c0(s)
        c_minus1 = self._c_minus1(s)
        c_minus2 = self._c_minus2(s)
        d1 = self._d1(s)
        return c0 * (c0**2 + d1 * c_minus1) / (c_minus1**2 - c0 * c_minus2) if c0 * (c0**2 + d1 * c_minus1) else 0.

    def _beta1(self, s: int) -> float:
        """
        Eq. C4
        :param s: Index of source layer
        """
        c0 = self._c0(s)
        c_minus1 = self._c_minus1(s)
        c_minus2 = self._c_minus2(s)
        d1 = self._d1(s)
        return (c0 * c_minus1 + d1 * c_minus2) / (c0 * c_minus2 - c_minus1**2) if (c0 * c_minus1 + d1 * c_minus2) else 0.

    def _beta2(self, s: int) -> float:
        """
        Eq. C5
        :param s: Index of source layer
        """
        c0 = self._c0(s)
        c_minus1 = self._c_minus1(s)
        c_minus2 = self._c_minus2(s)
        d1 = self._d1(s)
        return (c0**2 + d1 * c_minus1) / (c_minus1**2 - c0 * c_minus2) if (c0**2 + d1 * c_minus1) else 0.

    def _initial_estimate_q(self, s: int, X: float) -> float:
        """
        Eq. C1. Initial value of transformed ray parameter q for direct and
        reflected waves.
        :param s: Index of source layer
        :param X: Offset between source and receiver along x axis in m
        """
        alpha1 = self._alpha1(s)
        alpha2 = self._alpha2(s)
        beta1 = self._beta1(s)
        beta2 = self._beta2(s)
        numerator = (beta1 * X - alpha1 + sqrt((beta1**2 - 4 * beta2) * X**2
                                               + 2 * (2 * alpha2 - alpha1 * beta1) * X
                                               + alpha1**2))
        denominator = 2 * (alpha2 - beta2 * X)
        if denominator != 0:
            return numerator / denominator

    def _v_A(self) -> float:
        """
        Return v_A, the maximum velocity above the depth z(n-1)
        """
        velocities_bottom = self._velocity_model.velocities_bot[:-1]
        velocity_top = self._velocity_model.velocities_top[:-1]
        return max(np.max(velocities_bottom), np.max(velocity_top))

    def _v_M(self, source_position: np.ndarray, receiver_position: np.ndarray) -> float:
        """
        Return v_M, the highest velocity of the layers the ray passes through
        """
        # TODO This assumes direct ray between source and receiver so only the
        #  layers between them matter. This doesn't work for reflections and
        #  turning rays
        # first, get all the top and bottom velocities from the interfaces
        # between the source and receiver
        source_index = self.s
        receiver_index = self._velocity_model.layer_index(receiver_position[Index.Z])
        # sort indices because accessing array only works from low to high
        index_low = min(source_index, receiver_index)
        index_high = max(source_index, receiver_index)
        if index_low == index_high:
            # source and receiver are in the same layer
            return max(self._velocity_model.eval_at(*source_position),
                       self._velocity_model.eval_at(*receiver_position))
        velocities_bottom = self._velocity_model.velocities_bot[index_low:index_high]
        # +1 to exclude top velocity of the layer containing the upper point.
        # There are two cases:
        # - The point lays below the top of the layer, so the ray doesn't pass
        #   through the velocity at the top of the layer
        # - The point lays on the top depth, so the velocity should be included.
        #   This case in then handled by evaluating the velocity at both points.
        velocities_top = self._velocity_model.velocities_top[index_low+1:index_high+1]
        # second, also get velocity at the end points (source and receiver)
        v = max(self._velocity_model.eval_at(*source_position),
                self._velocity_model.eval_at(*receiver_position))
        return max(np.max(velocities_bottom), np.max(velocities_top), v)

    def _X_tilde(self, k: int, s: int, q: float) -> float:
        """
        Eq. A3
        :param k: Index of current layer
        :param s: Index of receiver layer
        :param q: Estimate of transformed ray parameter
        """
        return ((self._mu_tilde(k, s)
                * (sqrt(q**-2 + self._epsilon_tilde(k, s))
                   - sqrt(q**-2 + self._omega_tilde(k, s))) if self._delta_a(k) else 0.)
                + (self._h_tilde(k, s) / sqrt(q**-2 + self._epsilon_tilde(k, s)) if (1 - self._delta_a(k)) else 0.))

    def _X_tilde_prime(self, k: int, s: int, q: float) -> float:
        """
        Eq. B9
        :param k: Index of current layer
        :param s: Index of source layer
        :param q: Estimate of transformed ray parameter
        """
        return ((self._mu_tilde(k, s) / q**2
                * (1 / sqrt(1 + self._omega_tilde(k, s) * q**2)
                   - 1 / sqrt(1 + self._epsilon_tilde(k, s) * q**2)) if self._delta_a(k) else 0.)
                + (self._h_tilde(k, s) / (1 + self._epsilon_tilde(k, s) * q**2)**1.5 if (1 - self._delta_a(k)) else 0.))

    def _X_tilde_double_prime(self, k: int, s: int, q: float) -> float:
        """
        Eq. B10
        :param k: Index of current layer
        :param s: Index of source layer
        :param q: Estimate of transformed ray parameter
        """
        return ((self._mu_tilde(k, s) / q**3
                * ((2 + 3 * self._epsilon_tilde(k, s) * q**2)
                   / (1 + self._epsilon_tilde(k, s) * q**2)**1.5
                   - (2 + 3 * self._omega_tilde(k, s) * q**2)
                   / (1 + self._omega_tilde(k, s) * q**2)**1.5) if self._delta_a(k) else 0.)
                - (3 * self._h_tilde(k, s) * self._epsilon_tilde(k, s) * q
                / (1 + self._epsilon_tilde(k, s) * q**2)**2.5 if (1 - self._delta_a(k)) else 0.))

    def _f_tilde(self, s: int, q: float) -> float:
        """
        Eq. (16)
        :param s:
        :param q:
        :return:
        """
        return sum(self._X_tilde(k, s, q) for k in range(self.n+1)) - self.X

    def _f_tilde_prime(self, s: int, q: float) -> float:
        """
        Eq. B3
        :param s: Index of source layer
        :param q: Estimate for transformed ray parameter
        """
        return sum(self._X_tilde_prime(k, s, q) for k in range(self.n+1))

    def _f_tilde_double_prime(self, s: int, q: float) -> float:
        """
        Eq. B4
        :param s: Index of source layer
        :param q: Estimate of transformed ray parameter
        """
        return sum(self._X_tilde_double_prime(k, s, q) for k in range(self.n+1))

    def _next_q(self, s: int, q_old: float) -> float:
        A = 0.5 * self._f_tilde_double_prime(s, q_old)
        B = self._f_tilde_prime(s, q_old)
        C = self._f_tilde(s, q_old)
        q_plus = q_old + (-B + sqrt(B**2 - 4*A*C)) / (2*A)
        q_minus = q_old + (-B - sqrt(B**2 - 4*A*C)) / (2*A)
        f_tilde_plus = self._f_tilde(s, q_plus)
        f_tilde_minus = self._f_tilde(s, q_minus)
        if abs(f_tilde_plus) < abs(f_tilde_minus):
            return q_plus
        else:
            return q_minus

    @staticmethod
    def q_to_p(q: float, v_M: float) -> float:
        """
        Backtransform modified ray parameter q to standard ray parameter p.
        Reordered from eq. 15.
        :param q: modified ray parameter
        """
        return sqrt(q**2 / (v_M**2 + v_M**2 * q**2))

    def trace(self, source_position, receiver_position, accuracy: float = 0.0001) -> float:
        source_index = self._velocity_model.layer_index(source_position[Index.Z])
        # insert a and b of source layer as first item into a and b
        # this is done to fix inconsistencies of indexing, where the paper sums
        # over k=0...n while a only has n elements.
        # TODO array is mutable, so every call to this function mutates a and b
        #  make a copy
        self.a = np.insert(self.a, 0, self.a[source_index])
        self.b = np.insert(self.b, 0, self.b[source_index])
        # keep notation short and close to paper
        self.zs = source_position[Index.Z]
        self.s = source_index
        # maximum velocity above the depth z(n-1)
        self.v_A = self._v_A()
        # highest velocity of the layers the ray passes through
        self.v_M = self._v_M(source_position, receiver_position)
        self.X = abs(source_position[Index.X] - receiver_position[Index.X])

        q = self._initial_estimate_q(source_index, self.X)
        while True:
            q_next = self._next_q(source_index, q)
            if abs(q - q_next) < accuracy:
                q = q_next
                break
            q = q_next
        return self.q_to_p(q, self.v_M)


if __name__ == '__main__':
    vm = VelocityModel3D.from_file("/home/darius/git/double-beam/fang2019model.txt")
    twopoint = TwoPointRayTracing(vm)
    a = time.time()
    p = twopoint.trace((434., 0., 500.), (868., 0., 0))
    b = time.time()
    print(f"p: {p:.6e}")
    # compare to known good result
    assert p == 0.00024143895460092447, "Wrong result in raytracing"
    print("Runtime: ", b-a)
