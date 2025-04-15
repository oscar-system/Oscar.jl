kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
a, b, c, d = gens(kQ)

# compute Hartshorne's example from Section 20.5 in 24HLC (24 hours of local cohohomology)
# M = k[Q] (as a k[Q]-module)
I_M = ideal(kQ, [])
M = quotient_ring_as_module(I_M)
inj_res = injective_res(M, 3)

I = ideal(kQ, [a, b])

# cohomological degree 0
H0 = zeroth_local_cohomology(quotient_ring_as_module(I_M), I)
is_zero(H0)

# cohomological degree 1
H1 = local_cohomology(I_M, I, 1)
[h for h in H1.sectors if dim(h.H)>0] #sectors with non-zero local cohomomology
lc_zero(H1)

# cohomological degree 2 
H2 = local_cohomology(I_M, I, 2)
[h for h in H2.sectors if dim(h.H)>0] #sectors with non-zero local cohomology
lc_zero(H2)

#cohomological degree 3
H3 = local_cohomology(I_M, I, 3)
lc_zero(H3)
