
function lie_algebra_conformance_test(
  L::LieAlgebra{C}, parentT::DataType, elemT::DataType; num_random_tests::Int=10
) where {C<:RingElement}
  @testset "basic manipulation" begin
    x = L(rand(-10:10, dim(L)))

    @test parentT <: LieAlgebra{C}
    @test elemT <: LieAlgebraElem{C}
    @test L isa parentT
    @test x isa elemT

    @test parent_type(elemT) == parentT
    @test elem_type(parentT) == elemT

    @test parent(x) === L

    @test coefficient_ring(x) === coefficient_ring(L)
    @test elem_type(coefficient_ring(L)) == C

    # this block stays only as long as `ngens` and `gens` are not specialized for Lie algebras
    @test dim(L) == ngens(L)
    @test basis(L) == gens(L)
    @test all(i -> basis(L, i) == gen(L, i), 1:dim(L))

    @test dim(L) == length(basis(L))
    @test all(i -> basis(L, i) == basis(L)[i], 1:dim(L))

    @test iszero(zero(L))

    @test coefficients(x) == [coeff(x, i) for i in 1:dim(L)]
    @test all(i -> coeff(x, i) == x[i], 1:dim(L))
    @test sum(x[i] * basis(L, i) for i in 1:dim(L)) == x

    @test x == x
    @test deepcopy(x) == x
    @test hash(deepcopy(x)) == hash(x)
  end

  @testset "parent object call overload" begin
    @test L() == zero(L) == L(zeros(coefficient_ring(L), dim(L)))

    for _ in 1:num_random_tests
      coeffs = rand(-10:10, dim(L))
      x1 = L(coeffs)
      x2 = L(coefficient_ring(L).(coeffs))
      x3 = L(matrix(coefficient_ring(L), 1, dim(L), coeffs))
      x4 = L(sparse_row(matrix(coefficient_ring(L), 1, dim(L), coeffs)))
      x5 = L(x1)
      @test x1 == x2
      @test x1 == x3
      @test x1 == x4
      @test x1 == x5
    end
  end

  @testset "vector space axioms" begin
    for _ in 1:num_random_tests
      x = L(rand(-10:10, dim(L)))
      y = L(rand(-10:10, dim(L)))
      z = L(rand(-10:10, dim(L)))

      @test x + y == y + x
      @test x + (y + z) == (x + y) + z

      @test x + zero(L) == x
      @test zero(L) + x == x

      @test -x + x == zero(L)
      @test x + (-x) == zero(L)

      @test x - y == x + (-y)

      @test x * 0 == zero(L)
      @test 0 * x == zero(L)

      @test 2 * x == x + x
      @test x * 2 == x + x
      @test coefficient_ring(L)(2) * x == x + x
      @test x * coefficient_ring(L)(2) == x + x
    end
  end

  @testset "Lie algebra axioms" begin
    for _ in 1:num_random_tests
      x = L(rand(-10:10, dim(L)))
      y = L(rand(-10:10, dim(L)))
      z = L(rand(-10:10, dim(L)))

      @test x * y == bracket(x, y)

      @test (x + y) * z == x * z + y * z
      @test x * (y + z) == x * y + x * z

      @test x * x == zero(L)
      @test x * y == -(y * x)

      @test (x * (y * z)) + (y * (z * x)) + (z * (x * y)) == zero(L)
    end
  end
end

include("AbstractLieAlgebra-test.jl")
include("LinearLieAlgebra-test.jl")

@testset "LieAlgebras.LieAlgebra" begin
  # nothing here yet
end
