from numpy import empty, pi, meshgrid, linspace, sum, sin, exp
from numpy.fft import fftn, fftfreq, ifftn
f = open("conv.txt", "w")
for N in range(3, 50, 2):
    print "N =", N
    L = 2.
    x1d = linspace(0, L, N+1)[:-1]
    x, y, z = meshgrid(x1d, x1d, x1d)

    nr = 3*pi*exp(sin(pi*x)*sin(pi*y)*sin(pi*z))/4
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
f.close()
