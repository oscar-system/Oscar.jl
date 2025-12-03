@testset "solving" begin
    R, (x1,x2,x3) = polynomial_ring(QQ, [:x1, :x2, :x3])
    I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
    C, x = polynomial_ring(QQ, :x)
    elim = 84*x^4 - 40*x^3 + x^2 + x
    denom = 336*x^3 - 120*x^2 + 2*x + 1
    p1  = -184*x^3 + 80*x^2 - 4*x - 1
    p2 = -36*x^3 + 18*x^2 - 2*x
    prec = [32,48,64]
    for p in prec
        res = real_solutions(I, precision=p)
        sols = res[1]
        for g in gens(I)
            for s in sols
                @test abs(evaluate(g, s)) < QQFieldElem(2)^(-p+1)
            end
        end
    end
    @test res[2].vars == Symbol[:x1, :x2, :x3]
    @test res[2].elim == elim
    @test res[2].denom == denom
    @test res[2].param[1] == -1 * p1
    @test res[2].param[2] == -1 * p2

    rat_sols = Vector{QQFieldElem}[[1, 0, 0], [1//3, 0, 1//3]]
    rat_res = Oscar._rational_solutions(I)
    @test length(rat_res) == 2 && issetequal(rat_res, rat_sols)

    I = ideal(R,[x1-1,x1-1])
    @test_throws ErrorException real_solutions(I)
    @test_throws ErrorException Oscar._rational_solutions(I)
    # Issue #1040
    Qx, (x,) = polynomial_ring(QQ, [:x])
    I = ideal(Qx, [x^2 + 1])
    res = real_solutions(I)
    C, x = polynomial_ring(QQ, :x)
    @test res[1] == []
    @test res[2].vars == Symbol[:x]
    @test res[2].elim == x^2+1
    @test res[2].denom == 2*x

    # issue 1743
    R, (x1, x2) = polynomial_ring(QQ, [:x1, :x2])
    I = ideal(R, [x1 + ZZ(2)^100, x2 + ZZ(2)^100])
    sols = Vector{QQFieldElem}[[-1267650600228229401496703205376, -1267650600228229401496703205376]]
    @test sols == real_solutions(I)[1]
end

@testset "Rational solutions" begin
  R, (x, y) = QQ[:x, :y]
  I = ideal([x - 1, y - 1])
  J = ideal([x - 2, y - 3])
  pts = rational_solutions(I * J)
  @test length(pts) == 2
  @test issetequal(pts, Vector{QQFieldElem}[[1, 1], [2, 3]])

  k, a = quadratic_field(-1)
  R, (x, y) = k[:x, :y]
  I = ideal([x^2 + 1, y^3 - 1])
  pts = rational_solutions(I)
  @test length(pts) == 2
  @test issetequal(pts, Vector{elem_type(k)}[k.([a, 1]), k.([-a, 1])])

  k = GF(5)
  a = k(2)
  R, (x, y) = k[:x, :y]
  I = ideal([x^2 + 1, y^3 - 1])
  pts = rational_solutions(I)
  @test length(pts) == 2
  @test issetequal(pts, Vector{elem_type(k)}[k.([a, 1]), k.([-a, 1])])

  let # different coefficient ring
    k, a = quadratic_field(-1)
    R, (x, y) = QQ[:x, :y]
    I = ideal([x^2 + 1, y^3 - 1])
    pts = rational_solutions(I)
    @test length(pts) == 0
    pts = rational_solutions(k, I)
    @test issetequal(pts, Vector{elem_type(k)}[k.([a, 1]), k.([-a, 1])])
  end

  let # different coefficient ring
    k = GF(5)
    a = k(2)
    R, (x, y) = k[:x, :y]
    I = ideal([x^2 + 1, y^3 - 1])
    pts = rational_solutions(algebraic_closure(k), I)
    @test length(pts) == 6
  end
end

@testset "Rational solutions for homogeneous ideals" begin
  Q = projective_space(QQ, 2)
  x = homogeneous_coordinates(Q)
  i = ideal([x[1]-2*x[3], x[2]-3*x[3]])
  @test length(rational_solutions(i)) == 1

  let
    Q = projective_space(QQ, 2)
    x = homogeneous_coordinates(Q)
    i = ideal([x[1]^2-2*x[3]^2, x[2]-3*x[3]])
    @test length(rational_solutions(i)) == 0
    @test length(rational_solutions(algebraic_closure(QQ), i)) == 2
  end
end

@testset "Jacobian" begin
  R,(x,y) = polynomial_ring(algebraic_closure(QQ), [:x,:y])
  g = (x^2+y^2-1)*y*(x-3*y)
  si = ideal([g]) + jacobian_ideal(g)
  v=rational_solutions(si)
  @test length(v) == 5
end

@testset "Additional rational_solutions tests" begin
  let
    K = algebraic_closure(QQ)
    R, (x, y) = polynomial_ring(K, [:x, :y])

    I1 = ideal([x^2 - 1, y^2])
    I2 = ideal([x^2 - 1, (y - 2)^2])
    I = intersect(I1, I2)
    pts = rational_solutions(I)

    J1 = ideal([y^2 - 1, x^2])
    J2 = ideal([y^2 - 1, (x - 2)^2])
    I_swapped = intersect(J1, J2)
    pts_swapped = rational_solutions(I_swapped)

    @test issetequal(pts, [K.([1, 0]), K.([-1, 0]), K.([1, 2]), K.([-1, 2])])
    @test issetequal(pts_swapped, [K.([0, 1]), K.([0, -1]), K.([2, 1]), K.([2, -1])])
  end

  let
    K = algebraic_closure(QQ)
    R, (x, y) = polynomial_ring(K, [:x, :y])

    I = ideal([(x^2 - 1)^2, y^2 - 4])
    pts = rational_solutions(I)

    I_swapped = ideal([(y^2 - 1)^2, x^2 - 4])
    pts_swapped = rational_solutions(I_swapped)

    @test length(pts) == 4
    @test length(pts_swapped) == 4
  end

  let
    R, (x, y) = QQ[:x, :y]

    f1 = 2*x^2 - x*y + 2*y^2 - 2
    f2 = 2*x^2 - 3*x*y + 3*y^2 - 2
    I = ideal([f1, f2])
    pts = rational_solutions(I)

    f1_swapped = 2*y^2 - x*y + 2*x^2 - 2
    f2_swapped = 2*y^2 - 3*x*y + 3*x^2 - 2
    I_swapped = ideal([f1_swapped, f2_swapped])
    pts_swapped = rational_solutions(I_swapped)

    known = [(-1, 0), (-1//2, -1), (1//2, 1), (1, 0)]
    swapped_known = [(t[2], t[1]) for t in known]

    @test length(pts) == 4 &&
          all(t -> any(p -> p[1] == t[1] && p[2] == t[2], pts), known)

    @test length(pts_swapped) == 4 &&
          all(t -> any(p -> p[1] == t[1] && p[2] == t[2], pts_swapped), swapped_known)
  end

  let
    K = QQ
    R, (x, y) = polynomial_ring(K, [:x, :y])

    g = x*y*(y - 1)*(x - 1)
    I = ideal([g]) + jacobian_ideal(g)
    pts = rational_solutions(I)

    g_swapped = x*y*(x - 1)*(y - 1)
    I_swapped = ideal([g_swapped]) + jacobian_ideal(g_swapped)
    pts_swapped = rational_solutions(I_swapped)

    @test issetequal(pts, [K.([0, 1]), K.([1, 0]), K.([1, 1]), K.([0, 0])])
    @test issetequal(pts_swapped, [K.([1, 0]), K.([0, 1]), K.([1, 1]), K.([0, 0])])
  end

  let
    k = GF(101)
    R, (x, y) = polynomial_ring(k, [:x, :y])

    I1 = ideal([x^2, y^2 - 25])
    I2 = ideal([x^2 - 49, y^2])
    I = intersect(I1, I2)
    pts = rational_solutions(I)

    J1 = ideal([x^2 - 25, y^2])
    J2 = ideal([x^2, y^2 - 49])
    I_swapped = intersect(J1, J2)
    pts_swapped = rational_solutions(I_swapped)

    @test length(pts) == 4
    @test length(pts_swapped) == 4
  end

  let
    K = algebraic_closure(QQ)
    R, (u, v, w) = polynomial_ring(K, [:u, :v, :w])

    M = matrix(R, [
      0       u-1     v-1     w-1    u+v;
     -u+1       0    u+w   v-w    v;
     -v+1   -(u+w)    0    u-v    w;
     -w+1   -(v-w)  -(u-v)   0   u-w;
    -(u+v)   -v    -w   -(u-w)   0
    ])

    I = ideal(pfaffians(M, 4))
    pts = rational_solutions(I)

    @test length(pts) == degree(I)
  end

  let
    K = algebraic_closure(QQ)
    R, (x, y) = polynomial_ring(K, [:x, :y])

    f = x^3 + x*y + y^2 - 1
    g = x^2*y - y + 2
    I = ideal([f, g])
    pts = rational_solutions(I)

    I_swapped = ideal([y^3 + y*x + x^2 - 1, y^2*x - x + 2])
    pts_swapped = rational_solutions(I_swapped)

    @test length(pts) > 0 && length(pts) == degree(I)
    @test length(pts_swapped) > 0 && length(pts_swapped) == degree(I_swapped)
  end

  let
    k = GF(97)
    R, (x, y) = polynomial_ring(k, [:x, :y])

    I = ideal([x^3 + 2*y^2 + 3, y^3 + 5*x + 7])
    pts = rational_solutions(I)

    @test length(pts) > 0 && length(pts) == 1
  end
end
