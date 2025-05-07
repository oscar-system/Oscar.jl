@testset "first local cohomology example" begin
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)
x, y, z = gens(R_Q)

# get MonoidAlgebra
kQ = monoid_algebra(R_Q);
x, y, z = gens(kQ);

#example of computing local cohomology modules 
I_M = ideal(kQ, [x^2*z, x^4*y])

@test base_ring(I_M) == kQ

m = ideal(kQ, [x, y, z]) #maximal ideal
@test is_subset(I_M,m)

M = Oscar.quotient_ring_as_module(I_M)

#compute local cohomology
H0 = Oscar.zeroth_local_cohomology(quotient_ring_as_module(I_M),m)
@test !is_zero(H0)
H2 = 
H1 = Oscar.local_cohomology(I_M,m,1)
@test !Oscar.lc_zero(H1)
H1_sectors = [h for h in H1.sectors if dim(h.H) > 0]
@test all([dim(h.H) == 1 for h in H1_sectors])
@test all([ambient_dim(h.sector) == 2 for h in H1.sectors])

Oscar.local_cohomology(I_M,m,2)
@test Oscar.lc_zero(H2)

H3 = Oscar.local_cohomology(I_M,m,3)
@test Oscar.lc_zero(H3)

H4 = Oscar.local_cohomology(I_M,m,4)
@test Oscar.lc_zero(H4)
end

@testset "hartshorne example" begin
# this is Hartshorne's example from Section 20.5 in 24HLC (24 hours of local cohohomology)
kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
a, b, c, d = gens(kQ)

# M = k[Q] (as a k[Q]-module)
I_M = ideal(kQ, [])
M = Oscar.quotient_ring_as_module(I_M)
inj_res = Oscar.injective_res(M, 3)

I = ideal(kQ, [a, b])

# cohomological degree 0
H0 = Oscar.zeroth_local_cohomology(quotient_ring_as_module(I_M), I)
@test is_zero(H0)

# cohomological degree 1
H1 = Oscar.local_cohomology(I_M, I, 1)
H1_sectors = [h for h in H1.sectors if dim(h.H)>0] #sectors with non-zero local cohomomology
@test !Oscar.lc_zero(H1)
@test length(H1_sectors) == 1 

# cohomological degree 2 
H2 = Oscar.local_cohomology(I_M, I, 2)
H2_sectors = [h for h in H2.sectors if dim(h.H)>0] #sectors with non-zero local cohomology
@test !Oscar.lc_zero(H2)
@test length(H2_sectors) == 1

#cohomological degree 3
H3 = Oscar.local_cohomology(I_M, I, 3)
@test Oscar.lc_zero(H3)
end
