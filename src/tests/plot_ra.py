from pylab import plot, savefig, legend, xlim
from numpy import loadtxt, size

eigs = loadtxt("eigs.txt")
R = [float(x) for x in open("R.txt").read().split()]
#nmax = size(eigs, 1)
nmax = 5
for n in range(1, nmax+1):
    i = n-1
    plot(R, eigs[:, i], label="n = %d" % n)
legend()
xlim([0, 1])
savefig("eigs.png")
