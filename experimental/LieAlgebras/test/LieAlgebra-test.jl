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

  @testset "adjoint_matrix" begin
    @testset for n in 2:4, F in [QQ, GF(2), GF(3)]
      L = general_linear_lie_algebra(F, n)
      x = L(rand(-10:10, dim(L)))
      y = L(rand(-10:10, dim(L)))
      @test coefficients(x * y) == coefficients(y) * adjoint_matrix(x)
    end
  end

  @testset "Killing matrix" begin
    # Example from Hum72, Ch. 5.1
    F = QQ
    L = lie_algebra(
      F,
      2,
      [matrix(F, [0 1; 0 0]), matrix(F, [1 0; 0 -1]), matrix(F, [0 0; 1 0])],
      [:x, :h, :y],
    )
    @test killing_matrix(L) == matrix(F, [0 0 4; 0 8 0; 4 0 0])
  end

  @testset "Semisimplicity" begin
    L = special_linear_lie_algebra(QQ, 3)
    @test is_semisimple(L)
    L = special_orthogonal_lie_algebra(QQ, 3)
    @test is_semisimple(L)
    L = general_linear_lie_algebra(QQ, 3)
    @test !is_semisimple(L)

    L = special_linear_lie_algebra(GF(5), 3)
    @test is_semisimple(L)
    L = special_orthogonal_lie_algebra(GF(5), 3)
    @test is_semisimple(L)
    L = general_linear_lie_algebra(GF(5), 3)
    @test_broken !is_semisimple(L)
  end
end
