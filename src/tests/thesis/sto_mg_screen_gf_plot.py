from numpy import loadtxt, average, empty, array, shape, transpose, size
from pylab import semilogx, savefig, legend, title, grid, xlabel, ylabel, clf

data = loadtxt("mg_eigs_D.txt")
D = data[:, 0]
Etot = data[:, 1]

for i in range(6):
    semilogx(D, data[:, i+2], ":", label="E%d HF" % (i+1))
for i in range(6):
    semilogx(D, data[:, i+2+12], "-", label="E%d GF" % (i+1))
#for i in range(6):
#    semilogx(D, data[:, i+2], ":")
#semilogx(D, E1gf, "-", label="E1 GF")
#semilogx(D, E2gf, "-", label="E2 GF")
title("Dependence of eigenvalues on Debye length")
xlabel("D")
ylabel("Energy [a.u.]")
grid()
legend()
savefig("sto_mg_screen_gf.pdf")
