using Oscar

# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)

is_zero(modulus(R_Q))
x, y, z = gens(R_Q)

# get MonoidAlgebra
kQ = get_monoid_algebra(R_Q)
kQ = get_monoid_algebra_from_lattice([[0, 1], [1, 1], [2, 1]], QQ)
x, y, z = gens(kQ.algebra)

R = get_monoid_algebra_module(kQ, quotient_ring_as_module(ideal(R_Q, [])))

# define ideal over monoid algebra
I = ideal(kQ, [x^2*z, x^4*y])
J = ideal(kQ, [x^5*y, z^3])

M_I = quotient_ring_as_module(I.ideal)
M_J = quotient_ring_as_module(J.ideal)

_M = direct_sum(M_I, M_J; task=:none)
M = sub(_M, [y*_M[1]+y*_M[2], x^2*_M[2]])
N = get_monoid_algebra_module(kQ, M[1])

# irreducible resolution of M = kQ/I
# M = quotient_ring_as_module(I)
irr_res = irreducible_res(N)

# minimal injective resolution of kQ/I up to cohomological degree 3
inj_res = injective_res(N, 3)

inj_res.injMods[1].indecInjectives #J^0
inj_res.injMods[2].indecInjectives #J^1
inj_res.injMods[3].indecInjectives #J^2
inj_res.injMods[4].indecInjectives #J^3

inj_res.cochainMaps[1] #J^0 -> J^1
inj_res.cochainMaps[2] #J^1 -> J^2
inj_res.cochainMaps[3] #J^2 -> J^3

# irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_3 = inj_res.irrRes

# check if irreducible resolution
length(irr_res_3.irrSums)
image(irr_res_3.cochainMaps[1])[1] == kernel(irr_res_3.cochainMaps[2])[1]
image(irr_res_3.cochainMaps[2])[1] == kernel(irr_res_3.cochainMaps[3])[1]
image(irr_res_3.cochainMaps[3])[1] == kernel(irr_res_3.cochainMaps[4])[1]
image(irr_res_3.cochainMaps[4])[1] == kernel(irr_res_3.cochainMaps[5])[1]
image(irr_res_3.cochainMaps[5])[1] == kernel(irr_res_3.cochainMaps[6])[1]
is_injective(irr_res_3.inclusions[1])
is_injective(irr_res_3.inclusions[2])
is_injective(irr_res_3.inclusions[3])
is_injective(irr_res_3.inclusions[4])
is_injective(irr_res_3.inclusions[5])
is_injective(irr_res_3.inclusions[6])
is_surjective(irr_res_3.inclusions[6])
