# definition of polynomial ring k[x,y]
R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])

# get MonoidAlgebra
kQ = get_monoid_algebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ,[x^4,x^2*y^2,y^4])

# irreducible resolution of M = kQ/I
M = quotient_ring_as_module(I)
irr_res = irreducible_res(M)

# compute injective resolution up to cohomological degree 2
inj_res = injective_res(I,2)

inj_res.injMods[1].indecInjectives
inj_res.injMods[2].indecInjectives
inj_res.injMods[3].indecInjectives

inj_res.cochainMaps[1]
inj_res.cochainMaps[2]
inj_res.cochainMaps[3]


## irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_2 = inj_res.irrRes

#check if irreducible resolution
length(irr_res_2.irrSums)
image(irr_res_2.cochainMaps[1])[1] == kernel(irr_res_2.cochainMaps[2])[1]
image(irr_res_2.cochainMaps[2])[1] == kernel(irr_res_2.cochainMaps[3])[1]
is_injective(irr_res_2.inclusions[1])
is_injective(irr_res_2.inclusions[2])
is_injective(irr_res_2.inclusions[3])
is_surjective(irr_res_2.inclusions[3])