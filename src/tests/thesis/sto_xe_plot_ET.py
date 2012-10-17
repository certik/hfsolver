from numpy import loadtxt, average, empty, array, shape, transpose, size
from pylab import semilogy, savefig, legend, title, grid, xlabel, ylabel, plot

def get_differences(data):
    assert data.shape[1] == 2
    x = data[:, 0]
    y = data[:, 1]
    diffs = abs(y[1:] - y[:-1])
    return diffs

def guess_convergence(diffs, k=5):
    """
    Guess the converged values.

    "k" is the natural oscillation of the data. It needs to be higher
    than any artificial peaks in the convergence.
    """
    d = empty(len(diffs)-k)
    for i in range(len(diffs)-k):
        d[i] = average(diffs[i:i+k])
    i = 0
    while max(d[i+1:i+1+k]) < d[i]:
        i += 1
        if i+1+k > len(d): break
    return i + k

def get_errors(data, i):
    assert data.shape[1] == 2
    x = data[:, 0]
    y = data[:, 1]
    errs = abs(y - y[i])
    return errs

def convergence_plot(data, i):
    diffs = get_differences(data)
    errs = get_errors(data, i)
    x = data[:, 0]
    y = data[:, 1]
    semilogy(x[1:], diffs, "k--", label="STO abs(E-E_prev)")
    semilogy(x, errs, "k-", lw=2, label="STO abs(E-E_conv)")
    title("STO even-tempered basis convergence for Xe\nConverged value: dofs=%d, E_conv=%.11f" % (int(x[i]), y[i]))
    xlabel("DOFs")
    ylabel("Error [a.u.]")

#    plot(x, abs(yy - y[i]), "k^", label="STO even-tempered (ET)")
#    plot(x, abs(yy - y[i]), "k-")
    plot(x, abs(virial), "kx", label="STO virial theorem")
#    plot(x, abs(virial), "k--", label="STO ET virial theorem")
    #plot(x[1:], abs(yy[1:] - yy[:-1]), "k-.", label="STO ET error prev")
    #plot(x, abs(yy - yy[-1]), "k-", label="STO ET error")

    grid()
    legend()
    savefig("sto_xe_convergence_ET.pdf")


data = loadtxt("sto_xe_Etot.txt")
x = data[:, 1]
y = data[:, 2]
virial = data[:, 3]
data = array([x, y])
data = transpose(data)

diffs = get_differences(data)
i = guess_convergence(diffs, 5)
i = size(x)-1
print "Converged value: (%d, %.10f)" % tuple(data[i])
convergence_plot(data, i)
