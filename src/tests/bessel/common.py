from pylab import (plot, legend, savefig, clf, semilogy, grid, xlabel, ylabel,
    ylim, title)
from numpy import maximum, array

def make_plots(xx, yf, yrat, k):
    clf()
    xx = array(xx)
    yf = array(yf)
    yrat = array(yrat)
    title("Function values plot")
    plot(xx, yf, label="f(x)")
    plot(xx, yrat, label="rat")
    xlabel("x")
    ylabel("function value")
    legend()
    savefig("f%d.png" % k)

    clf()
    a = abs(yf-yrat)
    i = a.argmax()
    print "max abs error: %e" % a[i]
    print "x = %.16e" % xx[i]
    print "f(x): %.16e" % yf[i]
    print "r(x): %.16e" % yrat[i]
    b = abs(yf-yrat)/maximum(abs(yf), abs(yrat))
    i = b.argmax()
    print "max rel error: %e" % b[i]
    print "x = %.16e" % xx[i]
    print "f(x): %.16e" % yf[i]
    print "r(x): %.16e" % yrat[i]
    semilogy(xx, a, label="absolute error")
    semilogy(xx, b, label="relative error")
    grid()
    legend()
    title("Error plot")
    xlabel("x")
    ylabel("error")
    ylim([1e-18, 1e-14])
    savefig("error%d.png" % k)
