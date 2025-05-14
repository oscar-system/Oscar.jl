@testset "test is_normal" begin
kQ = monoid_algebra_from_lattice([[3,0,0,3],[2,1,0,3],[0,3,0,3],[3,0,1,0],[2,1,1,0],[0,3,1,0]],QQ)
@test !is_normal(kQ)
kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
@test is_normal(kQ)
kQ = monoid_algebra_from_lattice([[4,0],[3,1],[1,3],[0,4]],QQ)
@test !is_normal(kQ)
end
@testset "constuct MonoidAlgebras" begin
# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q, phi = quo(S, J)

# get MonoidAlgebra
kQ = monoid_algebra(R_Q)

# get same MonoidAlgebra from lattice
_kQ = monoid_algebra_from_lattice([[0, 1], [1, 1], [2, 1]], QQ)
f = hom(_kQ.algebra,kQ.algebra,[x,y,z]) #define isomorphism
@test is_injective(f) && is_surjective(f)
@test cone(kQ) == cone(_kQ)
@test dim(cone(kQ)) == 2


# definition of polynomial ring k[x,y]
R_Q, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]; weights=[[1, 0], [0, 1]])

# get MonoidAlgebra
kQ = monoid_algebra(R_Q)
_kQ = monoid_algebra_from_lattice([[1,0],[0,1]],QQ)
f = hom(_kQ.algebra,kQ.algebra,[x,y])
@test is_injective(f) && is_surjective(f)
@test dim(cone(kQ)) == 2
@test cone(kQ) == cone(_kQ)


#example with grading group ZZ^3
kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
a, b, c, d = gens(kQ)
@test dim(cone(kQ)) == 3
@test length([f for f in faces(kQ) if dim(f.poly) == 2]) == 4 && length(facets(cone(kQ))) == 4
@test is_pointed(kQ)
end

@testset "injective resolution over k[x,y]" begin
# definition of polynomial ring k[x,y]
R_Q, (x, y) = graded_polynomial_ring(QQ, ["x", "y"]; weights=[[1, 0], [0, 1]])

# get MonoidAlgebra
kQ = monoid_algebra(R_Q)

# define ideal over monoid algebra
I = ideal(kQ, [x^4, x^2*y^2, y^4])

#irreducible decomposition 
W = irreducible_decomposition(I)
@test I == intersect(W)

# irreducible resolution of M = kQ/I
M = quotient_ring_as_module(I)
irr_res = irreducible_resolution(M)
@test is_exact(irr_res.cochain_complex)

# compute injective resolution up to cohomological degree 2
inj_res = injective_resolution(I, 2)
@test inj_res.upto <= 2

## irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
inj_res_Q = inj_res.Q_graded_part
@test is_exact(inj_res_Q.cochain_complex)
end

@testset "injective resolutions over k[Q] = k[x,y,z]/(xz - y^2)" begin
# definition of monoid algebra as quotient of polynomial ring
S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
J = ideal(S, [x*z-y^2])
R_Q,_ = quo(S, J)

# get MonoidAlgebra
kQ = monoid_algebra(R_Q)



##first example
# define ideal over monoid algebra
I = ideal(kQ, [x^2*z, x^4*y])
J = ideal(kQ, [x^5*y, z^3])

M_I = quotient_ring_as_module(I)
M_J = quotient_ring_as_module(J)

_M = direct_sum(M_I, M_J; task=:none)
M,_ = sub(_M, [y*_M[1]+y*_M[2], x^2*_M[2]])

# compute irreducible resolution
irr_res = irreducible_resolution(M)
@test is_exact(irr_res.cochain_complex) #this test fails

# minimal injective resolution of kQ/I up to cohomological degree 3
inj_res = injective_resolution(M, 3)
@test inj_res.upto <= 3

# irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_Q = inj_res.Q_graded_part
@test is_exact(irr_res_Q.cochain_complex)



##second example
I = ideal(kQ, [y^3, x^3 * z])

#irreducible decomposition
W = irreducible_decomposition(I)
@test I == intersect(W...)

# compute irreducible resolution of M = kQ/I
M = quotient_ring_as_module(I)
irr_res = irreducible_resolution(M)

# compute injective resolution up to cohomological degree 3
inj_res = injective_resolution(I, 3)
@test inj_res.upto <= 3

# irreducible resolution that is the Q-graded part of the minimal injective resolution above
inj_res_Q = inj_res.Q_graded_part
@test is_exact(inj_res_Q.cochain_complex)



##third example
I = ideal(kQ, [x^2*z])

#irreducible decomposition
W = irreducible_decomposition(I)
@test I == intersect(W)

# compute irreducible resolution of M = kQ/I
M = quotient_ring_as_module(I)
irr_res = irreducible_resolution(M)
@test is_exact(irr_res.cochain_complex)

# compute injective resolution of kQ/I up to cohomological degree 3
inj_res = injective_resolution(I, 3)
@test inj_res.upto <= 3

# get irreducible resolution that is the Q-graded part of the minimal injective resolution above (shifted)
irr_res_Q = inj_res.Q_graded_part
@test is_exact(irr_res_Q.cochain_complex)



##fourth example 
I = ideal(kQ, [z^2, z*y, x^4])

#irreducible decomposition
W = irreducible_decomposition(I)
@test I == intersect(W...)

# irreducible resolution of M = R/I
M = quotient_ring_as_module(I)
irr_res = irreducible_resolution(M)
@test is_exact(irr_res.cochain_complex)

# injective resolution of M up to cohomological degree 3
inj_res = injective_resolution(I, 3)
@test inj_res.upto <= 3

irr_res_Q = inj_res.Q_graded_part
@test is_exact(irr_res_Q.cochain_complex)
end

@testset "injective resolution over ZZ^3-graded MonoidAlgebra" begin
kQ = monoid_algebra_from_lattice([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
a, b, c, d = gens(kQ)

# ### these examples take some time...
# #first example
# I = ideal(kQ, [a^2*b, c^2])

# W = irreducible_decomposition(I)
# @test I == intersect(W)

# inj_res = injective_resolution(I, 3)

# inj_res_Q = inj_res.Q_graded_part
# @test is_exact(inj_res_Q.cochain_complex)

# #second example
# I = ideal(kQ, [a^2*b, c^2, d*a^4])

# W = irreducible_decomposition(I)
# @test I == intersect(W)

# inj_res = injective_resolution(I, 3)
# @test is_exact(inj_res_Q.cochain_complex)
end
