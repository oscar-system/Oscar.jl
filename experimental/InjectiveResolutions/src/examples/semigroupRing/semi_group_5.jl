kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
R_Q = kQ.algebra
a, b, c, d = gens(kQ.algebra)
I = ideal(kQ, [a^2*b, c^2])
I = ideal(kQ, [a^2*b, c^2, d*a^4])
M = quotient_ring_as_module(I)
inj_res = injective_res(I, 3)

# compute Hartshorne's example from section 20.5 in 24HLC
Z = ideal(R_Q, [])
M = Oscar.MonoidAlgebra(kQ, quotient_ring_as_module(Z))
I_M = get_monoid_algebra_ideal(kQ, Z)
inj_res = injective_res(M, 3)
I = ideal(kQ, [a, b])

# cohomological degree 0
H0 = zeroth_local_cohomology(quotient_ring_as_module(I_M), I)

# cohomological degree 1
H1 = local_cohomology(I_M, I, 1)
lc_zero(H1)

# cohomological degree 2 
H2 = local_cohomology(I_M, I, 2)
lc_zero(H2)
