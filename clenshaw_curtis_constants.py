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

def clencurt(n1):
    """ Computes the Clenshaw Curtis nodes and weights """
    # http://www.scientificpython.net/pyblog/clenshaw-curtis-quadrature
    if n1 == 1:
        x = 0.
        w = 2.
    if n1 == 2:
        x = [-1., 1.]
        w = [1., 1.]
    else:
        n = n1 - 1
        C = np.zeros((n1,2))
        k = 2*(1+np.arange(np.floor(n/2)))
        C[::2,0] = 2/np.hstack((1, 1-k*k))
        C[1,1] = -n
        V = np.vstack((C,np.flipud(C[1:n,:])))
        F = np.real(np.fft.ifft(V, n=None, axis=0))
        x = F[0:n1,1]
        w = np.hstack((F[0,0],2*F[1:n,0],F[n,0]))
    return x, w

for i in range(1, 5):
    n = 2**i+ 1
    o_x = weights(n)[0][n /2 :]
    o_w = weights(n)[1][n /2 :]
    print "o_x", n, o_x
    print "o_w", n, o_w
    n = 2**(i+1)+ 1
    n_x = weights(n)[0][n /2 :]
    n_w = weights(n)[1][n /2 :]
    print "n_x", n, n_x
    print "n_w", n, n_w
    on_x = n_x[::2]
    on_w = n_w[::2]
    print "on_x", n, on_x
    print "on_w", n, on_w
    print "on_w/o_w", n, on_w/o_w
    # on_w/o_w is not a constant so we do need to hold the old f calls
    # looks like we can reuse x's but not w's
    print
