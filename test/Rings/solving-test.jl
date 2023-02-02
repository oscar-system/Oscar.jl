@testset "solving" begin
    R, (x1,x2,x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])
    I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
    C, x = PolynomialRing(QQ, "x")
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
                @test abs(evaluate(g, s)) < fmpq(2)^(-p+1)
            end
        end
    end
    @test res[2].vars == Symbol[:x1, :x2, :x3]
    @test res[2].elim == elim
    @test res[2].denom == denom
    @test res[2].param[1] == -1 * p1
    @test res[2].param[2] == -1 * p2

    rat_sols = Vector{fmpq}[[1, 0, 0], [1//3, 0, 1//3]]
    rat_res = Oscar._rational_solutions(I)
    @test length(rat_res) == 2 && issetequal(rat_res, rat_sols)

    I = ideal(R,[x1-1,x1-1])
    @test_throws ErrorException real_solutions(I)
    @test_throws ErrorException Oscar._rational_solutions(I)
    # Issue #1040
    Qx, (x,) = PolynomialRing(QQ, ["x"])
    I = ideal(Qx, [x^2 + 1])
    res = real_solutions(I)
    C, x = PolynomialRing(QQ, "x")
    @test res[1] == []
    @test res[2].vars == Symbol[:x]
    @test res[2].elim == x^2+1
    @test res[2].denom == 2*x

    # isssue 1743
    R, (x1, x2) = PolynomialRing(QQ, ["x1", "x2"])
    I = ideal(R, [x1 + fmpz(2)^100, x2 + fmpz(2)^100])
    sols = Vector{fmpq}[[-1267650600228229401496703205376, -1267650600228229401496703205376]]
    @test sols == real_solutions(I)[1]
end

@testset "Rational solutions" begin
  R, (x, y) = QQ["x", "y"]
  I = ideal([x - 1, y - 1])
  J = ideal([x - 2, y - 3])
  pts = rational_solutions(I * J)
  @test length(pts) == 2
  @test issetequal(pts, Vector{fmpq}[[1, 1], [2, 3]])

  k, a = quadratic_field(-1)
  R, (x, y) = k["x", "y"]
  I = ideal([x^2 + 1, y^3 - 1])
  pts = rational_solutions(I)
  @test length(pts) == 2
  @test issetequal(pts, Vector{elem_type(k)}[k.([a, 1]), k.([-a, 1])])

  k = GF(5)
  a = k(2)
  R, (x, y) = k["x", "y"]
  I = ideal([x^2 + 1, y^3 - 1])
  pts = rational_solutions(I)
  @test length(pts) == 2
  @test issetequal(pts, Vector{elem_type(k)}[k.([a, 1]), k.([-a, 1])])
end

@testset "Rational solutions for homogenous ideals" begin
  Q, x = proj_space(QQ, 2)
  i = ideal([x[1]-2*x[3], x[2]-3*x[3]])
  @test length(rational_solutions(i)) == 1
end
