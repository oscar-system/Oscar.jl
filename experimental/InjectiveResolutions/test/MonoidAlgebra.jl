@testset "Ring interface for MonoidAlgebra" begin
  a = 2
  function ConformanceTests.generate_element(A::MonoidAlgebra{<:FieldElem, T}) where {T<:MPolyQuoRing}
    return A(rand(base_ring(A.algebra), -a:a, 0:a, 0:a))
  end

  function ConformanceTests.generate_element(A::MonoidAlgebra{<:FieldElem, T}) where {T<:MPolyRing}
    return A(rand(A.algebra, -a:a, 0:a, 0:a))
  end

  S, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]; weights=[[0, 1], [1, 1], [2, 1]])
  J = ideal(S, [x*z-y^2])
  R_Q, phi = quo(S, J)

  # get MonoidAlgebra
  kQ = Oscar.MonoidAlgebra(R_Q)
  Oscar.ConformanceTests.test_Ring_interface(kQ)
  Oscar.ConformanceTests.test_Ring_interface_recursive(kQ)
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

