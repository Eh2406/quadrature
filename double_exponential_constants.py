# this is an experimental script to determine where the constants come from

import numpy as np


def weights(n):
    # Double Exponential Formulas for Numerical Integration By Hidetosi TAKAHASI* and Masatake MORI (3-3)
    # Integration(f(x), -1, 1) ~=
    # pi/2 * h * sum(f(tanh(pi/2 * sinh(n*h)))*(cosh(n*h)/ (cosh(pi/2 * sinh(n*h)))^2 ), n, -oo, oo)

    # so the ABCISSAS are:
    # tanh(pi/2 * sinh(n*h)

    # and the WEIGHTS are:
    # pi/2 * h *(cosh(n*h)/ (cosh(pi/2 * sinh(n*h)))^2

    # John D. Cook uses -oo = -3 and oo = 3

    u = np.linspace(-3, 3, n)
    h = u[1] - u[0]
    x = np.tanh(np.pi / 2 * np.sinh(u))
    w = np.pi / 2.0 * h * (np.cosh(u) / ((np.cosh(np.pi / 2 * np.sinh(u))) ** 2))
    return x, w


# and 1st layer has 7 fn cals so
n = 7
print weights(n)[0][n / 2:]
print weights(n)[1][n / 2:]

# and 2st layer has 6 fn cals more
n += 6
print weights(n)[0][n / 2:]
print weights(n)[1][n / 2:]

# and 3st layer has 12 fn cals more
n += 12
print weights(n)[0][n / 2:]
print weights(n)[1][n / 2:]


def test_with_sympy(n):
    import sympy
    a, b, c, x = sympy.symbols("a, b, c, x")
    f = a * x ** 2 + b * x + c
    print sum([w * f.subs(x, xi) for xi, w in zip(*weights(n))])
    print sympy.integrate(f, (x, -1, 1))


for i in range(1, 5):
    n = 2 ** i + 1
    o_x = weights(n)[0][n / 2:]
    o_w = weights(n)[1][n / 2:]
    print "o_x", n, o_x
    print "o_w", n, o_w
    n = 2 ** (i + 1) + 1
    n_x = weights(n)[0][n / 2:]
    n_w = weights(n)[1][n / 2:]
    print "n_x", n, n_x
    print "n_w", n, n_w
    on_x = n_x[::2]
    on_w = n_w[::2]
    print "on_x", n, on_x
    print "on_w", n, on_w
    print "on_w/o_w", n, on_w / o_w
    # on_w/o_w is a constant so we do not need to hold the old f calls
    print
