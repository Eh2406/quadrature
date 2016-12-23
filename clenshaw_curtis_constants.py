# http://www.sam.math.ethz.ch/~waldvoge/Papers/fejer.pdf
import numpy as np

# 2.4
def weights(n):
    n -= 1
    v = np.linspace(0, 1, n + 1) * np.pi
    j = np.array(range(1, n / 2 + 1))[:,None]

    b = 2. * np.ones(j.shape)
    b -= (j == n / 2.)

    c = 2 * np.ones(v.shape)
    c[0] -= 1
    c[-1] -= 1

    s = (b / (4. * j ** 2 - 1) * np.cos(2 * j * v)).sum(0)
    x = (np.cos(v) + 1.0) - 1.0
    w = c  / n * (1. - s)
    assert np.abs(w[0] - 1./(n ** 2. - 1. + (n % 2))) < 1e-12
    assert np.abs(w[-1] - 1./(n ** 2. - 1. + (n % 2))) < 1e-12
    return x, w

def test_with_sympy(n):
    import sympy;
    a, b, c, x = sympy.symbols("a, b, c, x")
    f = a* x **2 + b * x + c
    print sum([w*f.subs(x, xi) for xi, w in zip(*weights(n))])
    print sympy.integrate(f, (x, -1, 1))


for i in range(1, 5):
    n = 2**i+ 1
    print n, weights(n)[0][n /2 :]
    print n, weights(n)[1][n /2 :]
    print

# looks like we can reuse x's but not w
