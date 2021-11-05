@testset "msolve" begin
    R,(x,y) = PolynomialRing(QQ, ["x","y"])
    I = ideal(R,[x+13*y,x^2*y-12//5*y^2-x])
    singular_assure(I)
    res = Int32[2, 3], BigInt[1, 1, 13, 1, 1, 1, -12, 5, -1, 1], Int32[1, 0, 0, 1, 2, 1, 0, 2, 1, 0]
    @test Oscar.convert_singular_ideal_to_array(I.gens.S) == res
    R,(x,y) = PolynomialRing(FiniteField(32003), ["x","y"])
    I = ideal(R,[x+13*y,x^2*y-12*y^2-x])
    singular_assure(I)
    res = Int32[2, 3], Int32[1, 13, 1, 31991, 32002], Int32[1, 0, 0, 1, 2, 1, 0, 2, 1, 0]
    @test Oscar.convert_singular_ideal_to_array(I.gens.S) == res
    @test Singular.equal(I.gens.S, Oscar.convert_ff_gb_array_to_singular_ideal(Int32(2), res[1], res[2], res[3], base_ring(I.gens.S)))
    R, (x1,x2,x3,x4) = PolynomialRing(FiniteField(next_prime(2^28)), ["x1", "x2", "x3", "x4"])
    I = ideal(R,[x1+2*x2+2*x3+2*x4-1,
            x1^2+2*x2^2+2*x3^2+2*x4^2-x1,
            2*x1*x2+2*x2*x3+2*x3*x4-x2,
            x2^2+2*x1*x3+2*x2*x4-x3])

    f4(I);
    G = Singular.Ideal(I.gens.Sx,
            (I.gens.Sx).(
                         [x1 + 2*x2 + 2*x3 + 2*x4 - 1, x3^2 + 2*x2*x4 + 76695850*x3*x4 + 115043772*x4^2 + 115043768*x2 - 76695846*x3 - 38347924*x4, x2*x3 - 2*x2*x4 - 38347926*x3*x4 + 76695842*x4^2 - 57521884*x2 + 38347923*x3 - 115043767*x4, x2^2 + 2*x2*x4 - 115043767*x3*x4 - 38347921*x4^2 - 38347923*x2 + 115043768*x3 - 76695846*x4, x3*x4^2 - 29826161*x4^3 + 14913081*x2*x4 + 56338306*x3*x4 + 129246702*x4^2 - 4971027*x2 - 8285045*x3 - 9942054*x4, x2*x4^2 + 89478486*x4^3 + 29826162*x2*x4 - 4971027*x3*x4 - 29826162*x4^2 - 126761189*x2 + 9942054*x3, x4^4 + 35851649*x4^3 - 35550375*x2*x4 - 17256326*x3*x4 - 36956322*x4^2 - 76950494*x2 + 85955250*x3 - 101027337*x4]))
    @test Singular.equal(I.gb.S, G)

    R, (x1,x2,x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])
    I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
    C, x = PolynomialRing(QQ, "x")
    elim = 84*x^4 - 40*x^3 + x^2 + x
    denom = 336*x^3 - 120*x^2 + 2*x + 1
    p1  = -184*x^3 + 80*x^2 - 4*x - 1
    p2 = -36*x^3 + 18*x^2 - 2*x
    c1 = fmpz(-1)
    c2 = fmpz(-1)
    sols = Vector{fmpq}[[3197531698215911246794018079661//5070602400912917605986812821504, 1598765849107955623397009039829//5070602400912917605986812821504, -662230497759452443800611668909//5070602400912917605986812821504], [1, 0, 0], [143587366392252266220763399489//633825300114114700748351602688, 71793683196126133110381699745//633825300114114700748351602688, 173325283664805084153412401855//633825300114114700748351602688], [52818775009509558395695966891//158456325028528675187087900672, 3//2535301200456458802993406410752, 845100400152152934331135470251//2535301200456458802993406410752]]
    res = msolve(I)
    @test res[1][1] == elim
    @test res[1][2] == denom
    @test res[1][3][1] == p1
    @test res[1][3][2] == p2
    @test res[1][4][1] == c1
    @test res[1][4][2] == c2
    @test res[2] == sols

    I = ideal(R,[x1-1,x1-1])
    @test_throws ErrorException msolve(I)
    I = ideal(R,[x1-1,x1+1])
    @test msolve(I) == ([],[])
end
