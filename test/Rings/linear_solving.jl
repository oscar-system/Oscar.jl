@testset "linear solving" begin
  # MPolyRing
  let
    R, (x,y) = polynomial_ring(GF(3), [:x, :y])
    r() = rand(R, 1:3, 1:5, 1:2)
    n = rand(1:5)
    m = rand(1:5)
    l = rand(1:5)
    A = matrix(R, n, m, [r() for i in 1:n, j in 1:m])
    x = matrix(R, l, n, [r() for i in 1:l, j in 1:n])
    b = x * A
    @test @inferred can_solve(A, b)
    fl, y = @inferred can_solve_with_solution(A, b)
    @test fl && y * A == b
    fl, y, K = @inferred can_solve_with_solution_and_kernel(A, b)
    @test fl && y * A == b && is_zero(K * A)
  end

  # MPolyQuo
  let
    R, (x,y) = polynomial_ring(GF(3), [:x, :y])
    f = y^2+y+x^2
    C = ideal(R, [f])
    Q, = quo(R, C)
    r() = Q(rand(R, 1:2, 1:2, 1:2))
    n = rand(1:5)
    m = rand(1:5)
    l = rand(1:5)
    A = matrix(Q, n, m, [r() for i in 1:n, j in 1:m])
    x = matrix(Q, l, n, [r() for i in 1:l, j in 1:n])
    b = x * A
    @test @inferred can_solve(A, b)
    fl, y = @inferred can_solve_with_solution(A, b)
    @test fl && y * A == b
    fl, y, K = @inferred can_solve_with_solution_and_kernel(A, b)
    @test fl && y * A == b && is_zero(K * A)
  end

  # MPolyQuoLoc
  let
    R, (x, y, u, v) = GF(3)[:x, :y, :u, :v]
    f = x*v-y*u
    I = ideal(R, f)
    Q, p = quo(R, I)
    S = Oscar.MPolyComplementOfKPointIdeal(R, [QQ(1), QQ(0), QQ(1), QQ(0)])
    L, _ = localization(Q, S)
    a, b, c, d = L.((x, y, u, v)) 
    r() = L(Q(rand(R, 1:2, 1:2, 1:1)))
    n = rand(1:4)
    m = rand(1:4)
    l = rand(1:4)
    A = matrix(L, n, m, [r() for i in 1:n, j in 1:m])
    x = matrix(L, l, n, [r() for i in 1:l, j in 1:n])
    b = x * A
    @test @inferred can_solve(A, b)
    fl, y = @inferred can_solve_with_solution(A, b)
    @test fl && y * A == b
    fl, y, K = @inferred can_solve_with_solution_and_kernel(A, b)
    @test fl && y * A == b && is_zero(K * A)
  end

  # MPolyLoc
  let
    R, (x, y) = GF(7)[:x, :y]
    f = x^2+y^2-1
    I = ideal(R, f)
    S = Oscar.MPolyComplementOfPrimeIdeal(I)
    V, _ = localization(S)
    r() = V(rand(R, 1:2, 1:2, 1:1))
    n = rand(1:4)
    m = rand(1:4)
    l = rand(1:4)
    A = matrix(V, n, m, [r() for i in 1:n, j in 1:m])
    x = matrix(V, l, n, [r() for i in 1:l, j in 1:n])
    b = x * A
    @test @inferred can_solve(A, b)
    fl, y = @inferred can_solve_with_solution(A, b)
    @test fl && y * A == b
    fl, y, K = @inferred can_solve_with_solution_and_kernel(A, b)
    @test fl && y * A == b && is_zero(K * A)
  end

  # Laurent
  let
    F = GF(2)
    R, (x, y) = laurent_polynomial_ring(F, ["x", "y"])
    r() = rand(R, 1:2, 1:2, 1:1)
    n = rand(1:4)
    m = rand(1:4)
    l = rand(1:4)
    A = matrix(R, n, m, [r() for i in 1:n, j in 1:m])
    x = matrix(R, l, n, [r() for i in 1:l, j in 1:n])
    b = x * A
    @test @inferred can_solve(A, b)
    fl, y = @inferred can_solve_with_solution(A, b)
    @test fl && y * A == b
    fl, y, K = @inferred can_solve_with_solution_and_kernel(A, b)
    @test fl && y * A == b && is_zero(K * A)
  end
end
