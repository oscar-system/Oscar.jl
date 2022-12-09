@testset "solving" begin
    R, (x1,x2,x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])
    I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
    C, x = PolynomialRing(QQ, "x")
    elim = 84*x^4 - 40*x^3 + x^2 + x
    denom = 336*x^3 - 120*x^2 + 2*x + 1
    p1  = -184*x^3 + 80*x^2 - 4*x - 1
    p2 = -36*x^3 + 18*x^2 - 2*x
    sols = Vector{fmpq}[[744483363399261433351//1180591620717411303424, 372241681699630716673//1180591620717411303424, -154187553040555781639//1180591620717411303424], [1, 0, 0], [71793683196126133110381699745//316912650057057350374175801344, 71793683196126133110381699745//633825300114114700748351602688, 173325283664805084153412401855//633825300114114700748351602688], [196765270119568550571//590295810358705651712, 1//590295810358705651712, 196765270119568550571//590295810358705651712]]
    res = real_solutions(I)
    @test res[1] == sols
    @test res[2].vars == Symbol[:x1, :x2, :x3]
    @test res[2].elim == elim
    @test res[2].denom == denom
    @test res[2].param[1] == -1 * p1
    @test res[2].param[2] == -1 * p2

    rat_sols = Vector{fmpq}[[1, 0, 0], [1//3, 0, 1//3]]
    rat_res = Oscar._rational_solutions(I)
    @test rat_res == rat_sols

    I = ideal(R,[x1-1,x1-1])
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
