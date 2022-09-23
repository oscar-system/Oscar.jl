#R, x = PolynomialRing(QQ, ["x$i$j" for i in 1:3 for j in i:3])
#M = R[x[1] x[2] x[3]; x[2] x[4] x[5]; x[3] x[5] x[6]]
#E = QQ.([1, 0, 0, 1, 0, 1])
#G = AffineMatrixGroup(M, ideal(R, [x[2], x[3]]), E)

G = special_linear_group(2, QQ)
rho = canonical_representation(G)
rep = induced_representation_on_symmetric_power(G, 4)

O = omega_process(G)
I = nullcone_ideal(rep)
rop = reynolds_operator_from_omega_process(rep)
g = [rop(g) for g in gens(I)]

check = check_invariance_function(rep)
[check(h) for h in g]

A, m = invariant_ring(rep)

println("done.")
l = lift_to_invariant_polynomial_func(rep)
[l(4*y) for y in g]
[l(y^2) for y in g]
p = g[1]^3-5*g[2]^2 - 7*g[3]*g[1]
h = l(p)
h == A[1]^3-5*A[2]^2 - 7*A[3]*A[1]


#
#R, (z1, z2, z3, z4, x1, x2, x3) = QQ["z₁", "z₂", "z₃", "z₄", "x₁", "x₂", "x₃"]
#
#Z = zero(MatrixSpace(R, 3, 3))
#Z = R[z1^2 2*z1*z2 z2^2; z1*z3 z1*z4+z2*z3 z2*z4; z3^2 2*z3*z4 z4^2]
#X = R[x1; x2; x3]
#Y = Z*X
#S, x = QQ["x1", "x2", "x3"]
#p = x[2]^2 - x[1]*x[3]
#y = [Y[i, 1] for i in 1:3]
#P = evaluate(p, [x1, x2, x3])
#q = evaluate(p, y)
#D = Oscar.as_constant_differential_operator(det(Z))
