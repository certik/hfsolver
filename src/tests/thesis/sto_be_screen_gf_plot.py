from numpy import loadtxt, average, empty, array, shape, transpose, size
from pylab import semilogx, savefig, legend, title, grid, xlabel, ylabel, clf

data = loadtxt("eigs_D.txt")
D = data[:, 0]
Etot = data[:, 1]
E1 = data[:, 2]
E2 = data[:, 3]
E1gf = data[:, 4]
E2gf = data[:, 5]

semilogx(D, E1, "x", label="E1 HF")
semilogx(D, E2, "x", label="E2 HF")
semilogx(D, E1gf, "x", label="E1 GF")
semilogx(D, E2gf, "x", label="E2 GF")
title("Dependence of eigenvalues on Debye length")
xlabel("D")
ylabel("Energy [a.u.]")
grid()
legend()
savefig("sto_be_screen_gf.pdf")
