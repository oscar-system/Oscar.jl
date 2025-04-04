# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x * z - y^2])
R_Q, phi = quo(S, J)
x, y, z = gens(R_Q)

# get MonoidAlgebra
kQ = get_monoid_algebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ,[z^2, z*y, x^4])

# irreducible resolution of M = R/I
M = quotient_ring_as_module(I)
irr_res = irreducible_res(M)

# injective resolution of M up to cohomological degree 3
inj_res = injective_res(I,3)

inj_res.injMods[1].indecInjectives
inj_res.injMods[2].indecInjectives
inj_res.injMods[3].indecInjectives
inj_res.injMods[4].indecInjectives
inj_res.injMods[5].indecInjectives
inj_res.injMods[6].indecInjectives
inj_res.injMods[7].indecInjectives

inj_res.cochainMaps[1]
inj_res.cochainMaps[2]
inj_res.cochainMaps[3]
inj_res.cochainMaps[4]
inj_res.cochainMaps[5]
inj_res.cochainMaps[6]
inj_res.cochainMaps[7]

# irreducible resolution that is the Q-graded part of the minima injective resolution above (shifted)
irr_res_3 = inj_res.irrRes

# check if irreducible resolution
length(irr_res_3.irrSums)
image(irr_res_3.cochainMaps[1])[1] == kernel(irr_res_3.cochainMaps[2])[1]
image(irr_res_3.cochainMaps[2])[1] == kernel(irr_res_3.cochainMaps[3])[1] 
image(irr_res_3.cochainMaps[3])[1] == kernel(irr_res_3.cochainMaps[4])[1]
image(irr_res_3.cochainMaps[4])[1] == kernel(irr_res_3.cochainMaps[5])[1]
image(irr_res_3.cochainMaps[5])[1] == kernel(irr_res_3.cochainMaps[6])[1]
image(irr_res_3.cochainMaps[6])[1] == kernel(irr_res_3.cochainMaps[7])[1]
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

I = ideal(kQ,[z^3,z^2*y,x^4*z,x^3*y*z])
_M = quotient_ring_as_module(I)
M = shifted_module(I,[2,3])
