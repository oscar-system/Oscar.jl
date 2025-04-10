# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x * z - y^2])
R_Q, phi = quo(S, J)
x, y, z = gens(R_Q)

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ, [z^2, z*y, x^4])

# irreducible resolution of M = R/I
M = quotient_ring_as_module(I)
irr_res = irreducible_res(M)

# injective resolution of M up to cohomological degree 3
inj_res = injective_res(I, 3)

inj_res.inj_mods[1].indec_injectives
inj_res.inj_mods[2].indec_injectives
inj_res.inj_mods[3].indec_injectives
inj_res.inj_mods[4].indec_injectives
inj_res.inj_mods[5].indec_injectives
inj_res.inj_mods[6].indec_injectives
inj_res.inj_mods[7].indec_injectives

inj_res.cochain_maps[1]
inj_res.cochain_maps[2]
inj_res.cochain_maps[3]
inj_res.cochain_maps[4]
inj_res.cochain_maps[5]
inj_res.cochain_maps[6]
inj_res.cochain_maps[7]

# irreducible resolution that is the Q-graded part of the minima injective resolution above (shifted)
irr_res_3 = inj_res.irr_res

# check if irreducible resolution
length(irr_res_3.irr_sums)
image(irr_res_3.cochain_maps[1])[1] == kernel(irr_res_3.cochain_maps[2])[1]
image(irr_res_3.cochain_maps[2])[1] == kernel(irr_res_3.cochain_maps[3])[1]
image(irr_res_3.cochain_maps[3])[1] == kernel(irr_res_3.cochain_maps[4])[1]
image(irr_res_3.cochain_maps[4])[1] == kernel(irr_res_3.cochain_maps[5])[1]
image(irr_res_3.cochain_maps[5])[1] == kernel(irr_res_3.cochain_maps[6])[1]
image(irr_res_3.cochain_maps[6])[1] == kernel(irr_res_3.cochain_maps[7])[1]
is_injective(irr_res_3.inclusions[1])
is_injective(irr_res_3.inclusions[2])
is_injective(irr_res_3.inclusions[3])
is_injective(irr_res_3.inclusions[4])
is_injective(irr_res_3.inclusions[5])
is_injective(irr_res_3.inclusions[6])
is_injective(irr_res_3.inclusions[7])
is_injective(irr_res_3.inclusions[8])
is_injective(irr_res_3.inclusions[9])
is_surjective(irr_res_3.inclusions[9])

I = ideal(kQ, [z^3, z^2*y, x^4*z, x^3*y*z])
_M = quotient_ring_as_module(I)
M = shifted_module(I, [2, 3])
