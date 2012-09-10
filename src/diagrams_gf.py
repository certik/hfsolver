from numpy import empty, size

oo = 10000

def getn_L(L):
    n = -1
    for edge in L:
        if edge[0] != oo and edge[0] > n: n = edge[0]
        if edge[1] != oo and edge[1] > n: n = edge[1]
    return n

def L2Hugenholtz_index(L):
    def find_xy(n, L):
        f = []
        for edge in L:
            if edge[0] == n:
                f.append(edge[1])
        assert len(f) == 2
        x, y = f
        if x <= y:
            return x, y
        else:
            return y, x

    H = []
    for n in range(1, getn_L(L)+1):
        x, y = find_xy(n, L)
        H.append(x)
        H.append(y)
    return H

def Hugenholtz_index2L(H):
    def find_xy(H, i):
        f = []
        for n in range(len(H)):
            if H[n] == i:
                f.append(n)
        return f
    L = []
    for n in range(len(H)/2):
        L.append((n+1, H[2*n]))
        L.append((n+1, H[2*n+1]))
    for n in range(1, getn_L(L)+1):
        if len(find_xy(H, n)) == 1:
            if oo in H:
                L.append((0, n))
            else:
                L.append((oo, n))
    return L

def Goldstone_index2L(H, G):
    L = []
    for n in range(len(H)/2):
        L.append((n+1, 1, H[2*n], G[2*n]))
        L.append((n+1, 2, H[2*n+1], G[2*n+1]))
    return L

def Goldstone_index2Lindex(H, G, particles, holes):
    L = []
    pi = -1
    hi = -1
    for u, alpha, v, beta in Goldstone_index2L(H, G):
        if u > v:
            hi += 1
            index = holes[hi]
        else:
            pi += 1
            index = particles[pi]
        L.append((u, alpha, v, beta, index))
    return L, particles[:pi+1], holes[:hi+1]

def plus(H):
    n = size(H)/2
    return H + [n+1, n+1]

def c(H, j):
    L = Hugenholtz_index2L(H)
    n = getn_L(L)
    for ii, edge in enumerate(L):
        if edge[1] == edge[0] and edge[1] == n:
            i = ii+1
            break
    xi, yi = L[i-1]
    xj, yj = L[j-1]
    yi, yj = yj, yi
    L[i-1] = xi, yi
    L[j-1] = xj, yj
    return L2Hugenholtz_index(L)

def op1(H, j):
    return c(plus(H), j)

def op2(H, j, l):
    return c(op1(H, j), l)

def generate(Hlist):
    H2 = []
    for H in Hlist:
        for i in range(1, len(H)+4):
            new = op1(H, i)
            if new not in H2:
                H2.append(new)
            for j in range(1, len(H)+4):
                new = op2(H, i, j)
                if new not in H2:
                    H2.append(new)
    return H2

def bubbles(L):
    for edge in L:
        if edge[0] == edge[1]:
            return True
    return False

def filter_bubbles(Hlist):
    return [H for H in Hlist if not bubbles(Hugenholtz_index2L(H))]

def get_Goldstone_diagrams(H):
    n = len(H)/2
    def find_xy(H, i):
        f = []
        for n in range(len(H)):
            if H[n] == i:
                f.append(n)
        assert len(f) == 2
        return f
    def doit(G, i):
        if i == n+1:
            Gset.append(G)
            return
        x, y = find_xy(H, i)
        G1 = G[:]
        G1[x] = 1
        G1[y] = 2
        doit(G1, i+1)
        G2 = G[:]
        G2[x] = 2
        G2[y] = 1
        doit(G2, i+1)

    Gset = []
    doit([0]*len(H), 1)
    return Gset

def get_largest_v(edges):
    v2 = -1
    n2 = -1
    for n, (u, alpha, v, beta) in enumerate(edges):
        if v > v2:
            v2 = v
            n2 = n
    return n2

def get_loops(H, G):
    edges = Goldstone_index2L(H, G)
    loops = []
    while len(edges) > 0:
        loop = []
        n = get_largest_v(edges)
        u, alpha, v, beta = edges[n]
        loop.append(v)
        del edges[n]
        ok = True
        while ok:
            ok = False
            for n, (u2, alpha2, v2, beta2) in enumerate(edges):
                if v2 == u and beta2 == alpha:
                    del edges[n]
                    loop.append(v2)
                    u = u2
                    alpha = alpha2
                    ok = True
                    break
        loops.append(tuple(loop))
    loops.sort(key=lambda x: hash(x))
    return tuple(loops)

def filter_Goldstone_diagrams(Gset):
    # TODO: this works up to 3rd order, but get_loops() seem to do a few
    # mistakes in the 4rh order --- we need to investigate.
    G2set = []
    G2set_loops = []
    for H, G in Gset:
        loops = get_loops(H, G)
        if loops not in G2set_loops:
            G2set_loops.append(loops)
            G2set.append((H, G))
    return G2set

def n_holes(H):
    h = 0
    L = Hugenholtz_index2L(H)
    for edge in L:
        if edge[0] > edge[1]: h += 1
    return h

def expression(H, G, n_exchange):
    def find_xy(H, i):
        f = []
        for n in range(len(H)):
            if H[n] == i:
                f.append(n)
        assert len(f) == 2
        return f
    holes = ["a", "b", "c", "d", "e", "f", "g"]
    particles = ["r", "s", "t", "u", "v", "w", "x"]
    loops = get_loops(H, G)
    l = len(loops)
    h = n_holes(H)
    if n_exchange == 2:
        wg = 0.5
    else:
        wg = 1.0
    L, particles, holes = Goldstone_index2Lindex(H, G, particles, holes)
    int2 = []
    for n in range(len(H)/2):
        int2.append([0, 0, 0, 0])
    for u, alpha, v, beta, index in L:
        # v - in
        # u - out
        # alpha/beta: 1 (left), 2 (right)
        # <left in, right in | left out, right out>
        int2[v-1][beta-1] = index
        int2[u-1][alpha-1+2] = index
    en = []
    for n in range(len(H)/2-1):
        en.append([])
    for u, alpha, v, beta, index in L:
        for n in range(min(u,v)-1, max(u,v)-1):
            en[n].append(index)
    return holes, particles, l, h, wg, int2, en

def print_expression(holes, particles, l, h, wg, int2, en):
    print "(-1)^(%d+%d) * 2^%d" % (l, h, l),
    if wg == 0.5:
        print "* 1/2",
    for (i1, i2, i3, i4) in int2:
        print "<%s%s|%s%s>" % (i1, i2, i3, i4),
    for e in en:
        print "/ (",
        for i in e:
            if i in particles:
                print "-",
            else:
                print "+",
            print "eps_%s" % i,
        print ")",
    print

def fortran(f, Gset, subroutine):
    f.write("""\
real(dp) function %s(moint2, lam, Noccupied) result(E2)
real(dp), intent(in) :: moint2(:), lam(:)
integer, intent(in) :: Noccupied
""" % subroutine)
    H, G = Gset[0]
    holes, particles, l, h, wg, int2, en = expression(H, G, len(Gset))
    f.write("integer :: %s\n" % (", ".join(holes + particles)))
    f.write("E2 = 0\n")
    indent = 0
    for h in holes:
        f.write(" "*indent + "do %s = 1, Noccupied\n" % h)
        indent += 4
    for p in particles:
        f.write(" "*indent + "do %s = Noccupied+1, size(lam)\n" % p)
        indent += 4
    f.write(" "*indent + "E2 = E2 + (  &\n")
    for H, G in Gset:
        holes, particles, l, h, wg, int2, en = expression(H, G, len(Gset))
        f.write(" "*indent + "        + (-1)**(%d+%d) * 2**%d" %(l, h, l))
        if wg == 0.5:
            f.write(" / 2._dp")
        f.write("  &\n")
        for (i1, i2, i3, i4) in int2:
            # Important note:
            # the ijkl2intindex() uses (i1 i3 | i2 i4) = <i1 i2 | i3 i4>
            f.write(" "*indent + \
                "      * moint2(ijkl2intindex(%s, %s, %s, %s))  &\n" % \
                (i1, i3, i2, i4))
    f.write(" "*indent + "    ) / (  &\n")
    for n, e in enumerate(en):
        f.write(" "*indent + "    (")
        for i in e:
            if i in particles:
                f.write("-")
            else:
                f.write("+")
            f.write("lam(%s)" % i)
        f.write(")")
        if n + 1 < len(en):
            f.write(" *  &")
        else:
            f.write("  &")
        f.write("\n")
    f.write(" "*indent + "    )\n")
    for p in particles:
        indent -= 4
        f.write(" "*indent + "end do\n")
    for h in holes:
        indent -= 4
        f.write(" "*indent + "end do\n")
    f.write("end function\n")

def plot_Goldstone(H, G, filename):
    from pylab import plot, clf, title, savefig, xlim, ylim, arrow
    from numpy import array, sqrt, sum, abs, linspace, sin, cos
    from math import atan2
    def arc(xx, yy):
        A = array([xx[0], yy[0]])
        B = array([xx[1], yy[1]])
        plot(A[0], A[1], "ko")
        plot(B[0], B[1], "ko")
        AB = B-A
        AB_len = sqrt(sum(AB**2))
        AB = AB / AB_len
        pAB = array([AB[1], -AB[0]])
        AB_mid = (A+B)/2.
        R = AB_mid + pAB * AB_len/1.5
        r = sqrt(sum((A-R)**2))
        P_arrow = R - pAB * r
        # Draw the arc from A to B centered at R
        phi1 = atan2((A-R)[1], (A-R)[0])
        phi2 = atan2((B-R)[1], (B-R)[0])
        n = 100 # number of divisions
        phi = linspace(phi1, phi2, n)
        xx = r*cos(phi)+R[0]
        yy = r*sin(phi)+R[1]
        plot(xx, yy, "k-", lw=2)
        x, x2 = xx[n/2-1:n/2+1]
        y, y2 = yy[n/2-1:n/2+1]
        dx = x2-x
        dy = y2-y
        arrow(x, y, dx/2, dy/2, fc="black", head_width=0.05)
    loops = get_loops(H, G)
    L = Goldstone_index2L(H, G)
    clf()
    title("H = %r\nG = %r   loops = %r" % (H, G, loops))
    n = len(H)/2
    for i in range(1, n+1):
        plot([i, i], [1, 2], "k--")
    for u, alpha, v, beta in L:
        arc([u, v], [alpha, beta])
    xlim([0, n+1])
    ylim([0, 3])
    savefig(filename)

def plot_Hugenholtz(H, filename):
    from pylab import plot, clf, title, savefig, xlim, ylim, arrow
    from numpy import array, sqrt, sum, abs, linspace, sin, cos
    from math import atan2
    def arc(xx, yy, index, draw_dots=True):
        A = array([xx[0], yy[0]])
        B = array([xx[1], yy[1]])
        if draw_dots:
            plot(A[0], A[1], "ko")
            plot(B[0], B[1], "ko")
        AB = B-A
        AB_len = sqrt(sum(AB**2))
        AB = AB / AB_len
        pAB = array([AB[1], -AB[0]])
        AB_mid = (A+B)/2.
        if index == 1:
            R = AB_mid + pAB * AB_len/1.5
        else:
            R = AB_mid + pAB * AB_len/4.
        r = sqrt(sum((A-R)**2))
        P_arrow = R - pAB * r
        # Draw the arc from A to B centered at R
        phi1 = atan2((A-R)[1], (A-R)[0])
        phi2 = atan2((B-R)[1], (B-R)[0])
        n = 100 # number of divisions
        phi = linspace(phi1, phi2, n)
        xx = r*cos(phi)+R[0]
        yy = r*sin(phi)+R[1]
        plot(xx, yy, "k-", lw=2)
        x, x2 = xx[n/2-1:n/2+1]
        y, y2 = yy[n/2-1:n/2+1]
        dx = x2-x
        dy = y2-y
        arrow(x, y, dx/2, dy/2, fc="black", head_width=0.05)
    L = Hugenholtz_index2L(H)
    clf()
    title("H = %r" % H)
    n = len(H)/2
    cache = []
    for u, v in L:
        if (u, v) in cache:
            index = 2
        else:
            index = 1
            cache.append((u,v))
        if u == 0:
            if v == 1:
                plot([u, v], [1.5, 1.5], "k-", lw=2)
                arrow(0.5, 1.5, (v-u)/1000., 0, fc="black", head_width=0.05)
            else:
                arc([0, v], [1.5, 1.5], index, False)
        elif v == 0:
            if u == 1:
                plot([u, v], [1.5, 1.5], "k-", lw=2)
                arrow(0.5, 1.5, (v-u)/1000., 0, fc="black", head_width=0.05)
            else:
                arc([u, 0], [1.5, 1.5], index, False)
        elif u == oo or v == oo:
            if u == oo:
                if v == n:
                    plot([v, n+1], [1.5, 1.5], "k-", lw=2)
                    arrow(n+0.5, 1.5, -0.001, 0, fc="black", head_width=0.05)
                else:
                    arc([n+1, v], [1.5, 1.5], index, False)
            else:
                if u == n:
                    plot([u, n+1], [1.5, 1.5], "k-", lw=2)
                    arrow(n+0.5, 1.5, +0.001, 0, fc="black", head_width=0.05)
                else:
                    arc([u, n+1], [1.5, 1.5], index, False)
        else:
            arc([u, v], [1.5, 1.5], index)
    xlim([-1, n+2])
    ylim([0, 3])
    savefig(filename)


L1 = [
        (0, 1),
        (1, 2),
        (1, 3),
        (2, 1),
        (2, 3),
        (3, 2),
        (3, oo),
        ]

L2 = [
        (oo, 1),
        (1, 2),
        (1, 3),
        (2, 1),
        (2, 3),
        (3, 2),
        (3, 0),
        ]

Linit = [
        (0, 1),
        (1, 1),
        (1, oo),
        ]
Linitb = [
        (oo, 1),
        (1, 1),
        (1, 0),
        ]

H1 = [L2Hugenholtz_index(Linit), L2Hugenholtz_index(Linitb)]
H2 = generate(H1)
H3 = generate(H2)
#H4 = generate(H3)
H2 = filter_bubbles(H2)
H3 = filter_bubbles(H3)
#H4 = filter_bubbles(H4)
#H5 = generate(H4)

for n, H in enumerate(H1):
    print n, H, Hugenholtz_index2L(H)
print
for n, H in enumerate(H2):
    print n, H, Hugenholtz_index2L(H)
print
for n, H in enumerate(H3):
    print n, H, Hugenholtz_index2L(H)
    plot_Hugenholtz(H, "diag_gf%04d.png" % n)
