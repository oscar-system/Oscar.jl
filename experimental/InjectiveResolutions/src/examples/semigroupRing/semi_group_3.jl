# definition of monoid algebra as quotient of polynomial ring
S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
J = ideal(S,[x*z-y^2])
R_Q,phi = quo(S,J)
x,y,z = gens(R_Q)

# get MonoidAlgebra 
kQ = get_monoid_algebra(R_Q)

# define ideal over monoid algebra
I_Q = ideal(kQ,[x^2*z])

# compute irreducible resolution of M_Q = kQ/I_Q
M_Q = quotient_ring_as_module(I_Q)
irr_res = irreducible_res(M_Q)

# compute injective resolution of kQ/I_Q up to cohomological degree 3
inj_res = injective_res(I_Q,3)

inj_res.injMods[1].indecInjectives
inj_res.injMods[2].indecInjectives
inj_res.injMods[3].indecInjectives

inj_res.cochainMaps[1]
inj_res.cochainMaps[2]
inj_res.cochainMaps[3]

# get irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_3 = inj_res.irrRes

# check if irreducible resolution
length(irr_res_3.irrSums)
image(irr_res_3.cochainMaps[1])[1] == kernel(irr_res_3.cochainMaps[2])[1]
is_injective(irr_res_3.inclusions[1])
is_injective(irr_res_3.inclusions[2])
is_surjective(irr_res_3.inclusions[2])