# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)
x, y, z = gens(R_Q)

# get MonoidAlgebra 
kQ = Oscar.MonoidAlgebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ, [x^2*z])

#irreducible decomposition
W = irreducible_dec(I)
@test I == intersect(W...)

# compute irreducible resolution of M = kQ/I
M = quotient_ring_as_module(I)
irr_res = irreducible_res(M)
@test is_exact(irr_res.cochain_complex)

# compute injective resolution of kQ/I up to cohomological degree 3
inj_res = injective_res(I, 3)
@test inj_res.upto == 2

# get irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_Q = inj_res.Q_graded_part
@test is_exact(irr_res_Q.cochain_complex)

# check if irreducible resolution
@test image(irr_res_Q.cochain_maps[1])[1] == kernel(irr_res_Q.cochain_maps[2])[1]
@test is_injective(irr_res_Q.inclusions[1])
@test is_injective(irr_res_Q.inclusions[2])
@test is_injective(irr_res_Q.inclusions[3])
@test is_surjective(irr_res_Q.inclusions[3])
