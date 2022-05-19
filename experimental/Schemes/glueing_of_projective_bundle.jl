kk = GF(29)

# Set up the base ℙ¹ with coordinates s and t
R, (s,t) = PolynomialRing(kk, ["s", "t"])
S, _ = grade(R, [1, 1])

base_P1 = ProjectiveScheme(S)

# split this into the standard covering
base_covering = standard_covering(base_P1)

A1s = patches(base_covering)[1]
A1t = patches(base_covering)[2]

# Set up relative projective space of relative dimension 2 
# over both base patches
P2_s = projective_space(OO(A1s), ["xs", "ys", "zs"])

Cs = standard_covering(P2_s)

P2_t = projective_space(OO(A1t), ["xt", "yt", "zt"])

Ct = standard_covering(P2_t)

# Join the resulting schemes in a disjoint union with two 
# components
C = disjoint_union(Cs, Ct)

# Manually glue two (dense) patches of the two components
X = Cs[3]
Y = Ct[3]
x = gens(base_ring(OO(X)))
y = gens(base_ring(OO(Y)))
f = maximal_extension(X, Y, [x[1]//(x[3])^4, x[2]//(x[3])^6, 1//x[3]])
g = maximal_extension(Y, X, [y[1]//(y[3])^4, y[2]//(y[3])^6, 1//y[3]])
add_glueing!(C, Glueing(X, Y, restriction(f, domain(f), domain(g)), restriction(g, domain(g), domain(f))))

# Extend the glueing to the whole covered scheme
fill_transitions!(C)

X = CoveredScheme(C)
U = Cs[3]
(uw, vw, sw) = gens(base_ring(OO(U)))
fws = vw^2 - uw^3 - -14*sw^12 + 10*uw*sw^8 + 10*sw^5
I = IdealSheaf(X, [fws])
S = subscheme(I)

X = S
set_name!(X, "X")
H2 = IdealSheaf(X, [sw^2])
set_name!(H2, "H2")
D = WeilDivisor(H2)
L = IdealSheaf(X, [vw^2])
set_name!(L, "L")
D = D + 7*WeilDivisor(L)

I = IdealSheaf(X, [fws^8+1])
