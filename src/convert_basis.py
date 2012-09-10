from p631ss import basis_data

filename = "p631ss.txt"
f = open(filename, "w")
for Z in range(1, 92):
    if Z not in basis_data:
        continue
    f.write("%d\n" % Z)
    f.write("%d\n" % len(basis_data[Z]))
    for l, cgf in basis_data[Z]:
        f.write("%s %d\n" % (l, len(cgf)))
        for alpha, c in cgf:
            f.write("    %20.8f %20.8f\n" %(alpha, c))
print "Basis set saved into %s" % filename
