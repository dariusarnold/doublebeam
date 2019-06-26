"""
Two point ray tracing as described by
A fast and robust two-point ray tracing method in layered media with constant or
linearly varying layer velocity (Xinding Fang and Xiaofei Chen, 2019)
"""
from math import sqrt

import numpy as np

from doublebeam.core.raytracing.raytracing import VelocityModel1D


class MockVelocityModel:
    # This class mocks the velocity model shown in fig. 3 of Fang2019

    def __init__(self):
        # the velocity in every layer can be described as v_k = a_k*z + b_k
        # where k is the layer index starting at 0 for the top layer. a_k is the
        # velocity gradient and b_k the intercept.
        self.intercepts = np.array([1800., 2400., 2400., 2700., 2250.])
        self.gradients = np.array([4., 0., 1., 0., 1.5])
        self.depths = np.array([100., 200., 300., 400., 500])

    def layer_index(self, depth: float) -> int:
        """
        Return layer index that contains the given depth point.
        Layers are top inclusive, bottom exclusive, meaning a depth between two
        layers will give the index of the bottom layer.
        :param depth: Depth in m
        >>> vm = MockVelocityModel()
        >>> vm.layer_index(50.)
        0
        >>> vm.layer_index(499.9)
        4
        >>> vm.layer_index(300.)
        2
        """
        return np.searchsorted(self.depths, depth)

    def eval_at(self, depth: float) -> float:
        """
        Return velocity in m/s at the given depth
        :param depth: m
        >>> vm = MockVelocityModel()
        >>> vm.eval_at(0.)
        1800.0
        >>> vm.eval_at(500.)
        3000.0
        >>> vm.eval_at(350.)
        2700.0
        >>> vm.eval_at(50.)
        2000.0
        >>> vm.eval_at(100.)
        2200.0
        """
        index = self.layer_index(depth)
        return self.intercepts[index] + self.gradients[index] * depth



class TwoPointRayTracing:

    def __init__(self, velocity_model: MockVelocityModel):
        self._velocity_model: MockVelocityModel = velocity_model
        self.a = velocity_model.gradients
        self.b = velocity_model.intercepts
        self.z = velocity_model.depths
        self.n = len(velocity_model.depths)

    def _epsilon(self, k: int, s: int) -> float:
        """
        Eq. A5
        :param k: Index of current layer
        :param s: Index of source layer
        """
        if k == 0:
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
        return self._mu(k, s) * self.v_M / self.a[k] if self.a[k] != 0 else 0.

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
            return (self._delta_a(k) * self._mu_tilde(k, s)
                    * (sqrt(self._epsilon_tilde(k, s) - self._q_A() ** -2)
                       - sqrt(self._omega_tilde(k, s) - self._q_A() ** -2))
                    + (1 - self._delta_a(k)) * self._h_tilde(k, s)
                    / sqrt(self._epsilon_tilde(k, s) - self._q_A() ** -2))

        return (sum(d0_iter(k) for k in range(0, self.n-1))
                + self._mu_tilde(self.n, s) * sqrt(self._epsilon_tilde(self.n, s)
                                                   - self._q_A() ** -2))

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
            return (self._delta_a(k) * 0.5 * self._mu_tilde(k, s)
                    * (self._epsilon_tilde(k, s) - self._omega_tilde(k, s))
                    + (1 - self._delta_a(k)) * self._h_tilde(k, s))

        return sum(d1_iter(k) for k in range(0, self.n))

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
            return (self._delta_a(k) * self._mu_tilde(k, s)
                    * (self._delta_epsilon(k, s) * sqrt(self._epsilon_tilde(k, s))
                       - self._delta_omega(k, s) * sqrt(self._omega_tilde(k, s)))
                    + (1 - self._delta_a(k)) * self._delta_epsilon(k, s)
                    * self._h_tilde(k, s) / sqrt(self._epsilon_tilde(k, s)))

        return sum(c0_iter(k) for k in range(0, self.n))

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
            return ((1 - self._delta_a(k)) * (1 - self._delta_epsilon(k, s))
                    * self._h_tilde(k, s))

        return sum(c1_iter(k) for k in range(0, self.n))

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
            return (self._delta_a(k) * (self._delta_omega(k, s) - self._delta_epsilon(k, s))
                    * self._mu_tilde(k, s))

        return sum(c_minus1_iter(k) for k in range(0, self.n))

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
            return (self._delta_a(k) * 0.5 * self._mu_tilde(k, s)
                    * (self._delta_epsilon(k, s) / sqrt(self._epsilon_tilde(k, s))
                       - ((1 / sqrt(self._omega_tilde(k, s))) if self._delta_omega(k, s) != 0 else 0.))
                    - (1 - self._delta_a(k)) * self._delta_epsilon(k, s)
                    * 0.5 * self._h_tilde(k, s) / self._epsilon_tilde(k, s)**1.5)

        return sum(c_minus2_iter(k) for k in range(0, self.n))

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
        return c0 * (c0**2 + d1 * c_minus1) / (c_minus1**2 - c0 * c_minus2)

    def _beta1(self, s: int) -> float:
        """
        Eq. C4
        :param s: Index of source layer
        """
        c0 = self._c0(s)
        c_minus1 = self._c_minus1(s)
        c_minus2 = self._c_minus2(s)
        d1 = self._d1(s)
        return (c0 * c_minus1 + d1 * c_minus2) / (c0 * c_minus2 - c_minus1**2)

    def _beta2(self, s: int) -> float:
        """
        Eq. C5
        :param s: Index of source layer
        """
        c0 = self._c0(s)
        c_minus1 = self._c_minus1(s)
        c_minus2 = self._c_minus2(s)
        d1 = self._d1(s)
        return (c0**2 + d1 * c_minus1) / (c_minus1**2 - c0 * c_minus2)

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
        return (beta1 * X - alpha1 + sqrt((beta1**2 - 4 * beta2) * X**2 + 2
                                          * (2 * alpha2 - alpha1 * beta1)
                                          * X + alpha1**2)
                / (2 * (alpha2 - beta2 * X)))

    def _v_A(self) -> float:
        """
        Return v_A, the maximum velocity above the depth z(n-1)
        """
        velocities_bottom = (self.z * self.a + self.b)[:-1]
        velocity_top = (np.insert(self.z[:-1], 0, 0.) * self.a + self.b)[:-1]
        return max(np.max(velocities_bottom), np.max(velocity_top))

    def _v_M(self, source_position: np.ndarray, receiver_position: np.ndarray) -> float:
        """
        Return v_M, the highest velocity of the layers the ray passes through
        """
        # TODO This assumes direct ray between source and receiver so only the
        # layers between them matter. This doesn't work for reflections and
        # turning rays
        source_index = self._velocity_model.layer_index(source_position[2])
        receiver_index = self._velocity_model.layer_index(receiver_position[2])
        index_low, index_high = min(source_index, receiver_index), max(source_index, receiver_index)
        velocities_bottom = (self.a * self.z + self.b)[index_low:index_high]
        velocities_top = (np.insert(self.z[:-1], 0, 0.) *self.a + self.b)[index_low:index_high]
        v = max(self._velocity_model.eval_at(source_position[2]),
                self._velocity_model.eval_at(receiver_position[2]))
        return max(np.max(velocities_bottom), np.max(velocities_top), v)

    def _X_tilde(self, k: int, s: int, q: float) -> float:
        """
        Eq. A3
        :param k: Index of current layer
        :param s: Index of receiver layer
        :param q: Estimate of transformed ray parameter
        """
        return (self._delta_a(k) * self._mu_tilde(k, s)
                * (sqrt(q**-2 +  self._epsilon_tilde(k, s))
                   - sqrt(q**-2 + self._omega_tilde(k, s)))
                + ( 1 - self._delta_a(k)) * self._h_tilde(k, s)
                / sqrt(q**-2 + self._epsilon_tilde(k, s)))

    def _X_tilde_prime(self, k: int, s: int, q:float) -> float:
        """
        Eq. B9
        :param k: Index of current layer
        :param s: Index of source layer
        :param q: Estimate of transformed ray parameter
        """
        return (self._delta_a(k) * self._mu_tilde(k, s) / q**2
                * (1 / sqrt(1 + self._omega_tilde(k, s) * q**2)
                   - 1 / sqrt(1 + self._epsilon_tilde(k, s) * q**2))
                + (1 - self._delta_a(k)) * self._h_tilde(k, s)
                / (1 + self._epsilon_tilde(k, s) * q**2)**1.5)

    def _X_tilde_double_prime(self, k: int, s: int, q: float) -> float:
        """
        Eq. B10
        :param k: Index of current layer
        :param s: Index of source layer
        :param q: Estimate of transformed ray parameter
        """
        return (self._delta_a(k) * self._mu_tilde(k, s) / q**3
                * ((2 + 3* self._epsilon_tilde(k, s) * q**2)
                   / (1 + self._epsilon_tilde(k, s) * q**2)**1.5
                   - (2 + 3 * self._omega_tilde(k, s) * q**2)
                   / (1 + self._omega_tilde(k, s) * q**2)**1.5)
                - ( 1 - self._delta_a(k))
                * 3 * self._h_tilde(k, s) * self._epsilon_tilde(k, s) * q
                / (1 + self._epsilon_tilde(k, s) * q**2)**2.5)

    def _f_tilde(self, s: int, q: float) -> float:
        """
        Eq. (16)
        :param s:
        :param q:
        :return:
        """
        return sum(self._X_tilde(k, s, q) for k in range(self.n)) - self.X

    def _f_tilde_prime(self, s: int, q: float) -> float:
        """
        Eq. B3
        :param s: Index of source layer
        :param q: Estimate for transformed ray parameter
        """
        return sum(self._X_tilde_prime(k, s, q) for k in range(self.n))

    def _f_tilde_double_prime(self, s: int, q: float) -> float:
        """
        Eq. B4
        :param s: Index of source layer
        :param q: Estimate of transformed ray parameter
        """
        return sum(self._X_tilde_double_prime(k, s, q) for k in range(self.n))

    def _next_q(self, s: int, q_old: float) -> float:
        A = 0.5 * self._f_tilde_double_prime(s, q_old)
        B = self._f_tilde_prime(s, q_old)
        C = self._f_tilde(s, q_old)
        delta_q_plus = q_old + (-B + sqrt(B**2 -4*A*C)) / (2*A)
        delta_q_minus = q_old + (-B - sqrt(B**2 -4*A*C)) / (2*A)
        f_tilde_plus = self._f_tilde(s, delta_q_plus)
        f_tilde_minus = self._f_tilde(s, delta_q_minus)
        if abs(f_tilde_plus) < abs(f_tilde_minus):
            return delta_q_plus
        else:
            return delta_q_minus

    def trace(self, source_position, receiver_position):
        source_index = self._velocity_model.layer_index(source_position[2])
        # keep notation short and close to paper
        self.zs = source_position[2]
        s = source_index
        # maximum velocity above the depth z(n-1)
        self.v_A = self._v_A()
        # TODO Debug assert, remove this if successfully tested
        assert self.v_A == 2700, "Wrong value for v_A"
        # highest velocity of the layers the ray passes through
        self.v_M = self._v_M(source_position, receiver_position)

        self.X = abs(source_position[0] - receiver_position[0])

        q = self._initial_estimate_q(s, self.X)
        while True:
            print(q)
            q = self._next_q(s, q)



if __name__ == '__main__':
    vm = MockVelocityModel()
    twopoint = TwoPointRayTracing(vm)
    twopoint.trace((434., 0., 500.), (0., 868., 0))
