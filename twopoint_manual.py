import cmath
import math

import numpy as np
from math import sqrt, sin, radians

n = 5
a = np.array((4, 0, 1, 0, 1.5))
b = np.array((1800, 2400, 2400, 2700, 2250))
z = np.array((0, 100, 200, 300, 400, 500))

accuracy = 1E-10

source = (434., 0., 500.)
receiver = (868., 0., 0)
#source, receiver = receiver, source

zs = source[2]
s = np.searchsorted(z, zs, side="right")
if s > 5:
    s = 5

a = np.insert(a, 0, a[s-1])
b = np.insert(b, 0, b[s-1])

epsilon = np.array([(a[k] * z[k - 1] + b[k])**2 for k in range(1, n + 1)])
epsilon = np.insert(epsilon, 0, (a[s] * z[s - 1] + b[0])**2)

omega = np.array([(a[k] * z[k] + b[k])**2 for k in range(1, n + 1)])
omega = np.insert(omega, 0, (a[s] * zs + b[s])**2)

h = np.array([(a[k] * z[k-1] + b[k]) * (z[k]-z[k-1]) for k in range(1, n+1)])
h = np.insert(h, 0, (a[s] * z[s-1] + b[s]) * (zs - z[s-1]))


def mu_k(k, s):
    if 1 <= k <= s - 1:
        return 1
    if s <= k <= n:
        return 0
    if k == 0:
        return mu_k(s, s)

mu = np.array([mu_k(k, s) for k in range(0, n + 1)])
vM = 3000
vA = 2700

mu_tilde = np.divide(mu * vM, a, out=np.zeros_like(a), where=a!=0)

h_tilde = mu * h / vM

epsilon_tilde = 1 - epsilon / vM**2

omega_tilde = 1 - omega / vM**2

delta_a = np.abs(np.sign(a))

d1 = np.sum(delta_a * 0.5 * mu_tilde * (epsilon_tilde - omega_tilde) + (1 - delta_a) * h_tilde)

delta_epsilon = np.abs(np.sign(epsilon_tilde))

delta_omega = np.abs(np.sign(omega_tilde))

c0 = np.divide(delta_a * mu_tilde * (delta_epsilon * np.sqrt(epsilon_tilde)
                                     - delta_omega * np.sqrt(omega_tilde))
               + (1 - delta_a) * delta_epsilon * h_tilde, np.sqrt(epsilon_tilde), out=np.zeros_like(a), where=epsilon_tilde!=0)
c0 = np.sum(c0)

c1 = np.sum((1 - delta_a) * (1 - delta_epsilon) * h_tilde)

cminus1 = np.sum(delta_a * (delta_omega - delta_epsilon) * mu_tilde)

def safe_divide(a, b):
    """
    Divide a by b but return 0 on places where a is zero without dividing, since
    b could be zero there
    """
    return np.divide(a, b, out=np.zeros_like(a), where=b!=0)


cminus2 = (delta_a * 0.5 * mu_tilde * (safe_divide(delta_epsilon, np.sqrt(epsilon_tilde)) - safe_divide(delta_omega, np.sqrt(omega_tilde)))
           - (1 - delta_a) * delta_epsilon * h_tilde / (2 * epsilon_tilde**1.5))
cminus2 = np.sum(cminus2)

X = abs(source[0] - receiver[0])


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


def q_to_p(q: float, v_M: float) -> float:
    """
    Backtransform modified ray parameter q to standard ray parameter p.
    Reordered from eq. 15.
    :param q: modified ray parameter
    """
    return sqrt(q**2 / (v_M**2 + v_M**2 * q**2))


q = initial_q()
while True:
    print(q)
    q_next = next_q(q)
    if abs(q - q_next) < accuracy:
        q = q_next
        break
    q = q_next

p = q_to_p(q, vM)
p_expected = sin(radians(50)) / 3000
print("got: ", p, " expected: ", p_expected)
print(f"diff: {p - p_expected}")

assert sin(radians(50)) / 3000 == p, "fail"