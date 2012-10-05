from pylab import (plot, legend, savefig, clf, semilogy, grid, xlabel, ylabel,
    ylim, title)
from numpy import maximum, array

def make_plots(xx, yf, yrat):
    xx = array(xx)
    yf = array(yf)
    yrat = array(yrat)
    title("Function values plot")
    plot(xx, yf, label="f(x)")
    plot(xx, yrat, label="rat")
    xlabel("x")
    ylabel("function value")
    legend()
    savefig("f.png")

    clf()
    semilogy(xx, abs(yf-yrat), label="absolute error")
    semilogy(xx, abs(yf-yrat)/maximum(abs(yf), abs(yrat)),
            label="relative error")
    grid()
    legend()
    title("Error plot")
    xlabel("x")
    ylabel("error")
    ylim([1e-18, 1])
    savefig("error.png")
