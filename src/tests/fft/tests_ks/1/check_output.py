from json import load
from numpy import array

# Test results
D = load(open("output.json"))
eigs = array(D["md"][0]["energies"])
occ = array(D["md"][0]["occupation"])

# Reference results
Dref = load(open("output_ref.json"))
eigs_ref = array(Dref["md"][0]["energies"])
occ_ref = array(Dref["md"][0]["occupation"])

eig_err = abs(eigs-eigs_ref).max()
occ_err = abs(occ-occ_ref).max()
print("Eigs error:", eig_err)
print("Occupation error:", occ_err)
assert eig_err < 1e-6
assert occ_err < 1e-12
