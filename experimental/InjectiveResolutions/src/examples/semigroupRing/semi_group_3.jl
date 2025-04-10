# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)
x, y, z = gens(R_Q)

# get MonoidAlgebra 
kQ = Oscar.MonoidAlgebra(R_Q)

# define ideal over monoid algebra
I_Q = ideal(kQ, [x^2*z])

# compute irreducible resolution of M_Q = kQ/I_Q
M_Q = quotient_ring_as_module(I_Q)
irr_res = irreducible_res(M_Q)

# compute injective resolution of kQ/I_Q up to cohomological degree 3
inj_res = injective_res(I_Q, 3)

inj_res.inj_mods[1].indec_injectives
inj_res.inj_mods[2].indec_injectives
inj_res.inj_mods[3].indec_injectives

inj_res.cochain_maps[1]
inj_res.cochain_maps[2]
inj_res.cochain_maps[3]

# get irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_3 = inj_res.irr_res

# check if irreducible resolution
length(irr_res_3.irr_sums)
image(irr_res_3.cochain_maps[1])[1] == kernel(irr_res_3.cochain_maps[2])[1]
is_injective(irr_res_3.inclusions[1])
is_injective(irr_res_3.inclusions[2])
is_surjective(irr_res_3.inclusions[2])
