from numpy import empty, pi, meshgrid, linspace, sum
from numpy.fft import fftn, fftfreq, ifftn
E_exact = 128/(35*pi)
print "Hartree Energy (exact):      %.15f" % E_exact
f = open("conv.txt", "w")
for N in range(3, 384, 10):
    print "N =", N
    L = 2.
    x1d = linspace(0, L, N)
    x, y, z = meshgrid(x1d, x1d, x1d)

    nr = 3 * ((x-1)**2 + (y-1)**2 + (z-1)**2 - 1) / pi
    ng = fftn(nr) / N**3

    G1d = N * fftfreq(N) * 2*pi/L
    kx, ky, kz = meshgrid(G1d, G1d, G1d)
    G2 = kx**2+ky**2+kz**2
    G2[0, 0, 0] = 1  # omit the G=0 term

    tmp = 2*pi*abs(ng)**2 / G2
    tmp[0, 0, 0] = 0  # omit the G=0 term
    E = sum(tmp) * L**3
    print "Hartree Energy (calculated): %.15f" % E
    f.write("%d %.15f\n" % (N, E))

    Vh = 4*pi*ng / G2
    Vhr = ifftn(Vh).real
    print 0.5 * sum(Vhr*nr) * L**3
f.close()
