
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
  @testset "universal_enveloping_algebra" begin
    L = special_linear_lie_algebra(QQ, 2)

    U, L_to_U = universal_enveloping_algebra(L)
    e, f, h = basis(L)
    @test L_to_U(e) * L_to_U(f) - L_to_U(f) * L_to_U(e) == L_to_U(h)

    x = L(rand(-10:10, dim(L)))
    y = L(rand(-10:10, dim(L)))
    @test L_to_U(x) * L_to_U(y) - L_to_U(y) * L_to_U(x) == L_to_U(x * y)
  end

  @testset "Hum72, Exercise 2.2" begin
    @testset for n in 2:4, F in [QQ, GF(2)]
      L = general_linear_lie_algebra(F, n)

      derL = derived_algebra(L) # == sl_n

      @test dim(derL) == dim(L) - 1
      for b in basis(derL)
        @test tr(matrix_repr(b)) == 0
      end
    end
  end

  @testset "Hum72, Exercise 2.3" begin
    @testset for n in 2:4, F in [QQ, GF(2), GF(3)]
      L = general_linear_lie_algebra(F, n)
      cen = center(L) # == scalar matrices
      @test dim(cen) == 1
      b = matrix_repr(basis(cen, 1))
      @test divexact(b, b[1, 1]) == identity_matrix(F, n)

      L = special_linear_lie_algebra(F, n)
      cen = center(L)
      if is_divisible_by(n, characteristic(F))
        @test dim(cen) == 1
        b = matrix_repr(basis(cen, 1))
        @test divexact(b, b[1, 1]) == identity_matrix(F, n)
      else
        @test dim(cen) == 0
      end
    end
  end

  @testset "Hum72, Exercise 2.6" begin
    @testset for F in [QQ, GF(2), GF(3)]
      L = special_linear_lie_algebra(F, 3)

      @test_broken is_simple(L) == (characteristic(F) != 3)
    end
  end

  @testset "Perfectness" begin
    @testset for n in 2:5, F in [QQ, GF(2), GF(3)]
      L = general_linear_lie_algebra(F, n)
      @test !is_perfect(L)

      L = special_linear_lie_algebra(F, n)
      @test is_perfect(L) == !(n == 2 && characteristic(F) == 2)
    end
  end

  @testset "Solvability and nilpotency" begin
    L = abelian_lie_algebra(QQ, 3)
    @test is_nilpotent(L)
    @test is_solvable(L)

    @testset for n in 2:5, F in [QQ, GF(2), GF(3)]
      L = lie_algebra(
        F,
        n,
        [(b = zero_matrix(F, n, n); b[i, j] = 1; b) for i in 1:n for j in (i + 1):n],
        ["x_$(i)_$(j)" for i in 1:n for j in (i + 1):n],
      )
      @test is_nilpotent(L)
      @test is_solvable(L)
    end

    @testset for n in 2:5, F in [QQ, GF(2), GF(3)]
      L = lie_algebra(
        F,
        n,
        [(b = zero_matrix(F, n, n); b[i, j] = 1; b) for i in 1:n for j in i:n],
        ["x_$(i)_$(j)" for i in 1:n for j in i:n],
      )
      @test !is_nilpotent(L)
      @test is_solvable(L)
    end

    L = special_linear_lie_algebra(QQ, 3)
    @test !is_nilpotent(L)
    @test !is_solvable(L)
  end
end
