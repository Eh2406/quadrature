# this is an experimental script to determine where the constants come from

import numpy as np

# Double Exponential Formulas for Numerical Integration By Hidetosi TAKAHASI* and Masatake MORI (3-3)
# Integration(f(x), -1, 1) ~=
# pi/2 * h * sum(f(tanh(pi/2 * sinh(n*h)))*(cosh(n*h)/ (cosh(pi/2 * sinh(n*h)))^2 ), n, -oo, oo)

# so the ABCISSAS are:
# tanh(pi/2 * sinh(n*h)

# and the WEIGHTS are:
# pi/2 * h *(cosh(n*h)/ (cosh(pi/2 * sinh(n*h)))^2

# John D. Cook uses -oo = -3 and oo = 3

# and 1st layer has 7 fn cals so
u1 = np.linspace(-3, 3, 7)
h1 = u1[1] - u1[0]
ABCISSAS1 = np.tanh(np.pi/2 * np.sinh(u1))
print ABCISSAS1[int(len(ABCISSAS1)/2):]
WEIGHTS1 = np.pi/2.0 * h1 *(np.cosh(u1)/ ((np.cosh(np.pi/2 * np.sinh(u1)))**2))
print WEIGHTS1[int(len(WEIGHTS1)/2):]

# and 2st layer has 6 fn cals more
u2 = np.linspace(-3, 3, len(u1) + 6)
h2 = u2[1] - u2[0]
ABCISSAS2 = np.tanh(np.pi/2 * np.sinh(u2))
print ABCISSAS2[int(len(ABCISSAS2)/2):]
WEIGHTS2 = np.pi/2.0 * h2 *(np.cosh(u2)/ ((np.cosh(np.pi/2 * np.sinh(u2)))**2))
print WEIGHTS2[int(len(WEIGHTS2)/2):]

# and 3st layer has 12 fn cals more
u3 = np.linspace(-3, 3, len(u2) + 12)
h3 = u3[1] - u3[0]
ABCISSAS3 = np.tanh(np.pi/2 * np.sinh(u3))
print ABCISSAS3[int(len(ABCISSAS3)/2):]
WEIGHTS3 = np.pi/2.0 * h3 *(np.cosh(u3)/ ((np.cosh(np.pi/2 * np.sinh(u3)))**2))
print WEIGHTS3[int(len(WEIGHTS3)/2):]
