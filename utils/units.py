# SI units:
e = 1.602176565e-19     # C
kB = 1.3806488e-23      # J/K
eV = 1 * e                  # V * C = J/C * C = J
# For temperature, 1 eV/kB in K:
eV2K = eV / kB    # J / (J/K) = K

print "1e5 K =", 1e5 / eV2K, "eV"
print "15 keV =", 15e3 * eV2K / 1e6, "MK"
