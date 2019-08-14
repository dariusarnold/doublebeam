import math
from math import sqrt, sin, cos

import numpy as np

from doublebeam.models import VelocityModel3D
from doublebeam.utils import Index, safe_divide, horizontal_distance, angle_clockwise


class TwoPointRayTracing:
    """
    Implement two point ray tracing as described in
    A fast and robust two-point ray tracing method in layered media with
    constant or linearly varying layer velocity (Xinding Fang and Xiaofei Chen
    2019)
    This implementation only works for direct rays between the source and
    receiver. Turning rays and reflected rays are not implemented.
    """

    def __init__(self, velocity_model: VelocityModel3D):
        """
        Initialize two point ray tracing for the given velocity model
        """
        self._model = velocity_model
        self._a = velocity_model.gradients
        self._b = velocity_model.intercepts
        self._n = len(velocity_model)
        self._z = velocity_model.interface_depths

    def _v_M(self, source_position: np.ndarray,
             receiver_position: np.ndarray) -> float:
        """
        Return the highest velocity of the layers the direct ray passes through
        on the way from source to receiver.
        :param source_position: x, y, z triple of source coordinates
        :param receiver_position: x, y, z triple of receiver coordinates
        :return: Highest velocity (m/s) between source and receiver for a direct
        ray
        """
        receiver_index = self._model.layer_index(receiver_position)
        source_index = self._model.layer_index(source_position)
        if source_index == receiver_index:
            return max(self._model.eval_at(source_position),
                       self._model.eval_at(receiver_position))
        # first, get all the top and bottom velocities from the interfaces
        # between the source and receiver
        # sort indices because accessing array only works from low to high
        index_low = min(source_index, receiver_index)
        index_high = max(source_index, receiver_index)
        # dont +1 last index since the bottom of the bottom most layer is not
        # passed by the ray. If the ray starts there, evaluating the velocity
        # at source and receiver will catch it.
        velocities_bottom = self._model.velocities_bot[index_low:index_high]
        # +1 for first index to exclude top velocity of the layer containing the
        # upper point.
        # There are two cases:
        # - The point lays below the top of the layer, so the ray doesn't pass
        #   through the velocity at the top of the layer
        # - The point lays on the top depth, so the velocity should be included.
        #   This case in then handled by evaluating the velocity at both points.
        # +1 last index to include top velocity of the bottom most layer
        velocities_top = self._model.velocities_top[index_low+1:index_high+1]
        # second, also get velocity at the end points (source and receiver)
        v = max(self._model.eval_at(source_position),
                self._model.eval_at(receiver_position))
        return max(np.max(velocities_bottom), np.max(velocities_top), v)

    def trace(self, source: np.ndarray, receiver: np.ndarray,
              accuracy: float = 1E-10) -> np.ndarray:
        """
        Two point ray tracing from source to receiver position
        :param source: array of shape (3,) which contains the x y z coordinates
        of the source point
        :param receiver: array of shape (3,) which contains the x y z
        coordinates of the receiver point
        :param accuracy: Iteration for modified ray parameter q will be stopped
        if the difference between the current and the next value falls below
        this value.
        :return: Ray parameters px, py, pz in s/m
        """

        def initial_q():
            alpha1 = d1
            alpha2 = c0 * (c0**2 + d1 * cminus1) / (cminus1**2 - c0 * cminus2)
            beta1 = (c0 * cminus1 + d1 * cminus2) / (c0 * cminus2 - cminus1**2)
            beta2 = (c0**2 + d1 * cminus1) / (cminus1**2 - c0 * cminus2)
            numerator = (beta1 * X - alpha1 + sqrt((beta1**2 - 4 * beta2) * X**2
                                                   + 2 * (2 * alpha2 - alpha1 * beta1) * X
                                                   + alpha1**2))
            denominator = 2 * (alpha2 - beta2 * X)
            return numerator / denominator

        def X_tilde(q):
            return delta_a * mu_tilde * (np.sqrt(q**-2 + epsilon_tilde) - np.sqrt(q**-2 + omega_tilde)) \
                   + (1 - delta_a) * h_tilde / np.sqrt(q**-2 + epsilon_tilde)

        def X_tilde_prime(q):
            return delta_a * mu_tilde / (q**2) * (1 / np.sqrt(1 + omega_tilde * q**2)
                                                  - 1 / np.sqrt(1 + epsilon_tilde * q**2)) \
                   + (1 - delta_a) * h_tilde / (1 + epsilon_tilde * q**2)**1.5

        def X_tilde_double_prime(q):
            return delta_a * mu_tilde / q**3 * ((2 + 3 * epsilon_tilde * q**2) / (1 + epsilon_tilde * q**2)**1.5
                                                - (2 + 3 * omega_tilde * q**2) / (1 + omega_tilde * q**2)**1.5) \
                   - (1 - delta_a) * 3 * h_tilde * epsilon_tilde * q / (1 + epsilon_tilde * q**2)**2.5

        def f_tilde(q):
            return np.sum(X_tilde(q)) - X

        def f_tilde_prime(q):
            return np.sum(X_tilde_prime(q))

        def f_tilde_double_prime(q):
            return np.sum(X_tilde_double_prime(q))

        def next_q(q_old):
            A = 0.5 * f_tilde_double_prime(q)
            B = f_tilde_prime(q)
            C = f_tilde(q)
            delta_q_plus = (-B + sqrt(B**2 - 4 * A * C)) / (2 * A)
            delta_q_minus = (-B - sqrt(B**2 - 4 * A * C)) / (2 * A)
            q_plus = q_old + delta_q_plus
            q_minus = q_old + delta_q_minus
            if abs(f_tilde(q_plus)) < abs(f_tilde(q_minus)):
                return q_plus
            else:
                return q_minus

        top, bottom = self._model.vertical_boundaries()
        if not top <= source[Index.Z] <= bottom:
            raise ValueError(f"Source {source} outside of model")
        if not top <= receiver[Index.Z] <= bottom:
            raise ValueError(f"Receiver {receiver} outside of model")

        zs = source[Index.Z]
        a = self._a
        b = self._b
        n = self._n
        z = self._z

        # because paper uses 1 based indexing
        s = self._model.layer_index(zs) + 1
        # special case: bottom of the lowest layer should belong to the lowest
        # layer instead of the layer below as for the other layer
        if s > n:
            s = n

        a = np.insert(a, 0, a[s - 1])
        b = np.insert(b, 0, b[s - 1])

        # eq. A5
        epsilon = np.array([(a[k] * z[k - 1] + b[k])**2 for k in range(1, n + 1)])
        epsilon = np.insert(epsilon, 0, (a[s] * z[s - 1] + b[0])**2)

        # eq. A6
        omega = np.array([(a[k] * z[k] + b[k])**2 for k in range(1, n + 1)])
        omega = np.insert(omega, 0, (a[s] * zs + b[s])**2)

        # eq. A7
        h = np.array([(a[k] * z[k - 1] + b[k]) * (z[k] - z[k - 1]) for k in range(1, n + 1)])
        h = np.insert(h, 0, (a[s] * z[s - 1] + b[s]) * (zs - z[s - 1]))

        source_below_receiver = source[Index.Z] > receiver[Index.Z]
        mu = np.array([_mu_k(k, s, n, source_below_receiver) for k in range(0, n + 1)])
        vM = self._v_M(source, receiver)

        # eq. A10
        mu_tilde = np.divide(mu * vM, a, out=np.zeros_like(a), where=a != 0)

        # eq. A11
        h_tilde = mu * h / vM

        # eq. A12
        epsilon_tilde = 1 - epsilon / vM**2

        # eq. A13
        omega_tilde = 1 - omega / vM**2

        # eq.A9
        delta_a = np.abs(np.sign(a))

        # eq. C13
        d1 = np.sum(delta_a * 0.5 * mu_tilde * (epsilon_tilde - omega_tilde) + (1 - delta_a) * h_tilde)

        # eq. C18
        delta_epsilon = np.abs(np.sign(epsilon_tilde))

        # eq. C19
        delta_omega = np.abs(np.sign(omega_tilde))

        #eq. C14
        c0 = (delta_a * mu_tilde * (delta_epsilon * np.sqrt(epsilon_tilde) - delta_omega * np.sqrt(omega_tilde))
              + np.divide((1 - delta_a) * delta_epsilon * h_tilde, np.sqrt(epsilon_tilde), out=np.zeros_like(a),
                          where=epsilon_tilde != 0))
        c0 = np.sum(c0)

        # eq. C16
        cminus1 = np.sum(delta_a * (delta_omega - delta_epsilon) * mu_tilde)

        # eq. C17
        cminus2 = (delta_a * 0.5 * mu_tilde * (
                safe_divide(delta_epsilon, np.sqrt(epsilon_tilde)) - safe_divide(delta_omega, np.sqrt(omega_tilde)))
                   - safe_divide((1 - delta_a) * delta_epsilon * h_tilde, (2 * epsilon_tilde**1.5)))
        cminus2 = np.sum(cminus2)

        # horizontal distance between source and receiver
        X = horizontal_distance(source, receiver)
        q = initial_q()
        while True:
            q_next = next_q(q)
            if abs(q - q_next) < accuracy:
                q = q_next
                break
            q = q_next

        horizontal_slowness = _q_to_p(q, vM)
        c = self._model.eval_at(source)
        vertical_slowness = math.sqrt(c**-2 - horizontal_slowness**2)
        if source_below_receiver:
            # ray should travel upward from source in this case
            vertical_slowness *= -1
        receiver_surface_projection = receiver - source
        receiver_surface_projection[Index.Z] = 0
        x_axis = np.array((1, 0, 0))
        phi = angle_clockwise(x_axis, receiver_surface_projection)
        px = cos(phi) * horizontal_slowness
        py = sin(phi) * horizontal_slowness
        return np.array((px, py, vertical_slowness))


def _mu_k(k, s, n, source_below_receiver: bool = True) -> int:
    """
    Eq. A8
    :param k: Index of current layer
    :param s: Index of source layer
    :param n: Number of layers
    :param source_below_receiver: Fang2019s algorithm is only designed for
    sources below the receiver. If the source is above the receiver, this
    function has to be modified. Mu_k normally returns 1 for all layers above
    the source layer and 0 for all layers below. This must be reversed if the
    source is above the receiver
    """
    if 1 <= k <= s - 1:
        return 1 if source_below_receiver else 0
    if s <= k <= n:
        return 0 if source_below_receiver else 1
    if k == 0:
        return 1 - _mu_k(s, s, n, source_below_receiver)


def _q_to_p(q: float, v_M: float) -> float:
    """
    Backtransform modified ray parameter q to standard ray parameter p.
    Reordered from eq. 15.
    :param q: modified ray parameter
    """
    return sqrt(q**2 / (v_M**2 + v_M**2 * q**2))
