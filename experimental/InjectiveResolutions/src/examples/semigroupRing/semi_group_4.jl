# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x * z - y^2])
R_Q, phi = quo(S, J)
x, y, z = gens(R_Q)

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ, [z^2, z*y, x^4])

#irreducible decomposition
irreducible_dec(I)

# irreducible resolution of M = R/I
M = quotient_ring_as_module(I)
irr_res = irreducible_res(M)

# injective resolution of M up to cohomological degree 3
inj_res = injective_res(I, 3)
inj_res.upto

inj_res.inj_mods[1].indec_injectives #J^0
inj_res.inj_mods[2].indec_injectives #J^1
inj_res.inj_mods[3].indec_injectives #J^2
inj_res.inj_mods[4].indec_injectives #J^3

inj_res.cochain_maps[1] #J^0 -> J^1
inj_res.cochain_maps[2] #J^1 -> J^2
inj_res.cochain_maps[3] #J^2 -> J^3

# irreducible resolution that is the Q-graded part of the minima injective resolution above (shifted)
irr_res_Q = inj_res.irr_res

# check if irreducible resolution
length(irr_res_Q.irr_sums)
image(irr_res_Q.cochain_maps[1])[1] == kernel(irr_res_Q.cochain_maps[2])[1]
image(irr_res_Q.cochain_maps[2])[1] == kernel(irr_res_Q.cochain_maps[3])[1]
image(irr_res_Q.cochain_maps[3])[1] == kernel(irr_res_Q.cochain_maps[4])[1]
image(irr_res_Q.cochain_maps[4])[1] == kernel(irr_res_Q.cochain_maps[5])[1]
image(irr_res_Q.cochain_maps[5])[1] == kernel(irr_res_Q.cochain_maps[6])[1]
image(irr_res_Q.cochain_maps[6])[1] == kernel(irr_res_Q.cochain_maps[7])[1]
image(irr_res_Q.cochain_maps[7])[1] == kernel(irr_res_Q.cochain_maps[8])[1]
image(irr_res_Q.cochain_maps[8])[1] == kernel(irr_res_Q.cochain_maps[9])[1]
is_injective(irr_res_Q.inclusions[1])
is_injective(irr_res_Q.inclusions[2])
is_injective(irr_res_Q.inclusions[3])
is_injective(irr_res_Q.inclusions[4])
is_injective(irr_res_Q.inclusions[5])
is_injective(irr_res_Q.inclusions[6])
is_injective(irr_res_Q.inclusions[7])
is_injective(irr_res_Q.inclusions[8])
is_injective(irr_res_Q.inclusions[9])
is_surjective(irr_res_Q.inclusions[9])
