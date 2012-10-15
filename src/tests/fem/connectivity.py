# 2D connectivity table for tensor-product order p elements

from numpy import zeros

# Polynomial order of elements
p = 2
# Number of elements in x and y directions
nex = 2
ney = 3
# Boundary condition: 1 = Neumann, 2 = Dirichlet
ibc = 2

#================================================================

# Construct array of global nodes
# 0 = no associated basis function as in e.g. Dirichlet BCs
nodes = zeros((ney*p+1, nex*p+1))
inode = 0
for iy in range(ney*p+1):
    for ix in range(nex*p+1):
        if ibc == 2 and (ix == 0 or ix == nex*p or iy == 0 or iy == ney*p):
            continue
        inode += 1
        nodes[iy, ix] = inode
print "Global nodes:"
print nodes

# Construct connectivity table of global nodes in each element
gn = zeros((ney*nex, (p+1)**2))
iel = -1
for iey in range(ney):
    for iex in range(nex):
        iel += 1 # element number
        ix0 = iex*p # lower left corner
        iy0 = iey*p
        iln = -1
        # get global node numbers in element
        for iy in range(iy0, iy0+p+1):
            for ix in range(ix0, ix0+p+1):
                iln += 1 # local node number
                gn[iel, iln] = nodes[iy, ix] # global node number
print
print "Connectivity table (global node numbers in each element):"
print gn
