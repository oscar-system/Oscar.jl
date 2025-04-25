using Oscar

# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)
kQ = monoid_algebra_from_lattice([[0, 1], [1, 1], [2, 1]], QQ)
x, y, z = gens(kQ.algebra)

R = quotient_ring_as_module(ideal(kQ, []))

# define ideal over monoid algebra
I = ideal(kQ, [x^2*z, x^4*y])
J = ideal(kQ, [x^5*y, z^3])

M_I = quotient_ring_as_module(I)
M_J = quotient_ring_as_module(J)

_M = direct_sum(M_I, M_J; task=:none)
M_mod,_ = sub(_M, [y*_M[1]+y*_M[2], x^2*_M[2]])
M = M_mod
# M = monoid_algebra_module(kQ, M_mod[1])

# irreducible resolution of M = kQ/I
# M = quotient_ring_as_module(I)
# irr_res = irreducible_res(M)
irr_res = irreducible_res(M)

# minimal injective resolution of kQ/I up to cohomological degree 3
inj_res = injective_res(M, 3)

inj_res.inj_mods[1].indec_injectives #J^0
inj_res.inj_mods[2].indec_injectives #J^1
inj_res.inj_mods[3].indec_injectives #J^2
inj_res.inj_mods[4].indec_injectives #J^3

inj_res.cochain_maps[1] #J^0 -> J^1
inj_res.cochain_maps[2] #J^1 -> J^2
inj_res.cochain_maps[3] #J^2 -> J^3

# irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_3 = inj_res.irr_res

# check if irreducible resolution
length(irr_res_3.irr_sums)
image(irr_res_3.cochain_maps[1])[1] == kernel(irr_res_3.cochain_maps[2])[1]
image(irr_res_3.cochain_maps[2])[1] == kernel(irr_res_3.cochain_maps[3])[1]
image(irr_res_3.cochain_maps[3])[1] == kernel(irr_res_3.cochain_maps[4])[1]
image(irr_res_3.cochain_maps[4])[1] == kernel(irr_res_3.cochain_maps[5])[1]
image(irr_res_3.cochain_maps[5])[1] == kernel(irr_res_3.cochain_maps[6])[1]
image(irr_res_3.cochain_maps[6])[1] == kernel(irr_res_3.cochain_maps[7])[1]
image(irr_res_3.cochain_maps[7])[1] == kernel(irr_res_3.cochain_maps[8])[1]
image(irr_res_3.cochain_maps[8])[1] == kernel(irr_res_3.cochain_maps[9])[1]
image(irr_res_3.cochain_maps[9])[1] == kernel(irr_res_3.cochain_maps[10])[1]
is_injective(irr_res_3.inclusions[1])
is_injective(irr_res_3.inclusions[2])
is_injective(irr_res_3.inclusions[3])
is_injective(irr_res_3.inclusions[4])
is_injective(irr_res_3.inclusions[5])
is_injective(irr_res_3.inclusions[6])
is_injective(irr_res_3.inclusions[7])
is_injective(irr_res_3.inclusions[8])
is_injective(irr_res_3.inclusions[9])
is_injective(irr_res_3.inclusions[10])
is_surjective(irr_res_3.inclusions[10])
