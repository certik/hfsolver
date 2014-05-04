# This script provides a function rat(x) that gives rational approximation to
# the function besseli(4+0.5, x) / exp(x) for any x > 0. It plots an error plot
# showing, that both absolute and relative errors are less than 1e-15 against
# the exact answer calculated using mpmath.  The rational approximation was
# calculated using Mathematica, a sample script is given in the comment below.

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

from numpy import linspace
from math import exp
from sympy.mpmath import besseli, besselk

from common import make_plots
from hfsolver.special import Inu, Knu

xx = linspace(1e-10, 40, 10000)
for k in [0, 1, 2, 3, 4]:
    print "Doing:", k
    yf = [(besseli(k+0.5, x) / exp(x)) for x in xx]
    yrat = [Inu(k, x) for x in xx]
    make_plots(xx, yf, yrat, k)
