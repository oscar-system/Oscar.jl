kQ = get_monoid_algebra_from_lattice([[1,0,0],[1,1,0],[1,1,1],[1,0,1]],QQ)
R_Q = kQ.algebra
a,b,c,d = gens(kQ.algebra)
I = ideal(kQ,[a^2*b,c^2])
I = ideal(kQ,[a^2*b,c^2,d*a^4])
M = quotient_ring_as_module(I)
inj_res = injective_res(I,3)

Z = ideal(R_Q,[])
M = get_monoid_algebra_module(kQ,quotient_ring_as_module(Z))
inj_res = injective_res(M,3)