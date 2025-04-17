# definition of polynomial ring k[x,y]
R_Q, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]; weights=[[1, 0], [0, 1]])

# get MonoidAlgebra
kQ = Oscar.MonoidAlgebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ, [x^4, x^2*y^2, y^4])

#irreducible decomposition 
irr_dec = irreducible_dec(I)

# irreducible resolution of M = kQ/I
M = quotient_ring_as_module(I)
irr_res = irreducible_res(M)

# compute injective resolution up to cohomological degree 2
inj_res = injective_res(I, 2)

inj_res.inj_mods[1].indec_injectives
inj_res.inj_mods[2].indec_injectives
inj_res.inj_mods[3].indec_injectives

inj_res.cochain_maps[1]
inj_res.cochain_maps[2]
inj_res.cochain_maps[3]

## irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
inj_res_Q = inj_res.irr_res

#check if irreducible resolution
length(inj_res_Q.irr_sums)
image(inj_res_Q.cochain_maps[1])[1] == kernel(inj_res_Q.cochain_maps[2])[1]
image(inj_res_Q.cochain_maps[2])[1] == kernel(inj_res_Q.cochain_maps[3])[1]
is_injective(inj_res_Q.inclusions[1])
is_injective(inj_res_Q.inclusions[2])
is_injective(inj_res_Q.inclusions[3])
is_surjective(inj_res_Q.inclusions[3])
