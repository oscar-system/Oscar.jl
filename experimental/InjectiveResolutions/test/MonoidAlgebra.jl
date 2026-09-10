@testset "Ring interface for MonoidAlgebra" begin
  a = 2
  function ConformanceTests.generate_element(A::MonoidAlgebra{<:FieldElem, T}) where {T<:MPolyQuoRing}
    return A(rand(base_ring(A.algebra), -a:a, 0:a, 0:a))
  end

  function ConformanceTests.generate_element(A::MonoidAlgebra{<:FieldElem, T}) where {T<:MPolyRing}
    return A(rand(A.algebra, -a:a, 0:a, 0:a))
  end

  # get MonoidAlgebra
  kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)
  Oscar.ConformanceTests.test_Ring_interface(kQ)
  #Oscar.ConformanceTests.test_Ring_interface_recursive(kQ) # Needs AA update
end

@testset "constuct MonoidAlgebras" begin
  # get MonoidAlgebra
  kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]],QQ)

  kQ = monoid_algebra([[1,0],[0,1]],QQ)

  #example with grading group ZZ^3
  kQ = monoid_algebra([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
  a, b, c, d = gens(kQ)
  @test dim(cone(kQ)) == 3
  @test length([f for f in faces(kQ) if dim(f.poly) == 2]) == 4 && length(facets(cone(kQ))) == 4
  @test is_pointed(kQ)
end

@testset "minimal_generating_set for MonoidAlgebraIdeal" begin
  kQ = monoid_algebra([[0, 1], [1, 1], [2, 1]], QQ)
  x, y, z = gens(kQ)

  # x^2*z = x*y^2 in kQ, so x^2*z is redundant in (x*y, x^2*z)
  I = ideal(kQ, [x*y, x^2*z])
  mg = minimal_generating_set(I)
  @test all(g -> g in I, mg)
  @test ideal(kQ, mg) == I
  @test length(mg) < length(gens(I))

  # zero ideal
  I0 = ideal(kQ, elem_type(kQ)[])
  @test isempty(minimal_generating_set(I0))
end

@testset "test is_normal" begin
  kQ = monoid_algebra([[3,0,0,3],[2,1,0,3],[0,3,0,3],[3,0,1,0],[2,1,1,0],[0,3,1,0]],QQ)
  @test !is_normal(kQ)
  kQ = monoid_algebra([[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]], QQ)
  @test is_normal(kQ)
  kQ = monoid_algebra([[4,0],[3,1],[1,3],[0,4]],QQ)
  @test !is_normal(kQ)
end

