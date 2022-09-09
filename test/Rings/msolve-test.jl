@testset "msolve" begin
    R,(x,y) = PolynomialRing(QQ, ["x","y"])
    I = ideal(R,[x+13*y,x^2*y-12//5*y^2-x])
    singular_assure(I)
    res = Int32[2, 3], BigInt[1, 1, 13, 1, 1, 1, -12, 5, -1, 1], Int32[1, 0, 0, 1, 2, 1, 0, 2, 1, 0]
    @test Oscar.convert_singular_ideal_to_array(I.gens.S) == res
    R,(x,y) = PolynomialRing(GF(32003), ["x","y"], ordering=:degrevlex)
    I = ideal(R,[x+13*y,x^2*y-12*y^2-x])
    singular_assure(I)
    res = Int32[2, 3], Int32[1, 13, 1, 31991, 32002], Int32[1, 0, 0, 1, 2, 1, 0, 2, 1, 0]
    @test Oscar.convert_oscar_ideal_to_array(I) == res
    @test Oscar.convert_singular_ideal_to_array(I.gens.S) == res
    @test Singular.equal(I.gens.S, Oscar.convert_ff_gb_array_to_singular_ideal(Int32(2), res[1], res[2], res[3], base_ring(I.gens.S)))
    @test I.gens.O == Oscar.convert_ff_gb_array_to_oscar_array(Int32(2), res[1], res[2], res[3], base_ring(I), 0)
    R, (x1,x2,x3,x4) = PolynomialRing(GF(next_prime(2^28)), ["x1", "x2", "x3", "x4"], ordering=:degrevlex)
    I = ideal(R,[x1+2*x2+2*x3+2*x4-1,
            x1^2+2*x2^2+2*x3^2+2*x4^2-x1,
            2*x1*x2+2*x2*x3+2*x3*x4-x2,
            x2^2+2*x1*x3+2*x2*x4-x3])
    H = f4(I);
    G = gfp_mpoly[x1 + 2*x2 + 2*x3 + 2*x4 + 268435458
                  x3^2 + 2*x2*x4 + 76695850*x3*x4 + 115043772*x4^2 + 115043768*x2 + 191739613*x3 + 230087535*x4
                  x2*x3 + 268435457*x2*x4 + 230087533*x3*x4 + 76695842*x4^2 + 210913575*x2 + 38347923*x3 + 153391692*x4
                  x2^2 + 2*x2*x4 + 153391692*x3*x4 + 230087538*x4^2 + 230087536*x2 + 115043768*x3 + 191739613*x4
                  x3*x4^2 + 238609298*x4^3 + 14913081*x2*x4 + 56338306*x3*x4 + 129246702*x4^2 + 263464432*x2 + 260150414*x3 + 258493405*x4
                  x2*x4^2 + 89478486*x4^3 + 29826162*x2*x4 + 263464432*x3*x4 + 238609297*x4^2 + 141674270*x2 + 9942054*x3
                  x4^4 + 35851649*x4^3 + 232885084*x2*x4 + 251179133*x3*x4 + 231479137*x4^2 + 191484965*x2 + 85955250*x3 + 167408122*x4]
    @test H == G
    @test isdefined(I, :gb)
    @test I.gb[degrevlex(gens(base_ring(I)))].O == G
    H = f4(I, eliminate=2);
    G = gfp_mpoly[x3^2*x4 + 73209671*x3*x4^2 + 260301051*x4^3 + 188447115*x3^2 + 167207272*x3*x4 + 120660383*x4^2 + 210590781*x3 + 109814506*x4
                  x3^3 + 156877866*x3*x4^2 + 59264971*x4^3 + 224858274*x3^2 + 183605206*x3*x4 + 130731555*x4^2 + 110395535*x3 + 158620953*x4
                  x4^4 + 167618101*x3*x4^2 + 102789335*x4^3 + 193931678*x3^2 + 156155981*x3*x4 + 60823186*x4^2 + 239040667*x3 + 127377432*x4
                  x3*x4^3 + 99215126*x3*x4^2 + 261328123*x4^3 + 132228634*x3^2 + 93598185*x3*x4 + 85654356*x4^2 + 3613010*x3 + 240673711*x4]
    @test H == G
    @test I.gb[degrevlex(gens(base_ring(I))[3:end])].O == G

    R, (x1,x2,x3) = PolynomialRing(QQ, ["x1", "x2", "x3"])
    I = ideal(R, [x1+2*x2+2*x3-1, x1^2+2*x2^2+2*x3^2-x1, 2*x1*x2+2*x2*x3-x2])
    C, x = PolynomialRing(QQ, "x")
    elim = 84*x^4 - 40*x^3 + x^2 + x
    denom = 336*x^3 - 120*x^2 + 2*x + 1
    p1  = -184*x^3 + 80*x^2 - 4*x - 1
    p2 = -36*x^3 + 18*x^2 - 2*x
    c1 = fmpz(-1)
    c2 = fmpz(-1)
    sols = Vector{fmpq}[[744483363399261433351//1180591620717411303424, 372241681699630716673//1180591620717411303424, -154187553040555781639//1180591620717411303424], [1, 0, 0], [71793683196126133110381699745//316912650057057350374175801344, 71793683196126133110381699745//633825300114114700748351602688, 173325283664805084153412401855//633825300114114700748351602688], [196765270119568550571//590295810358705651712, 1//590295810358705651712, 196765270119568550571//590295810358705651712]]
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
    # Issue #1040
    Qx, (x,) = PolynomialRing(QQ, ["x"])
    I = ideal(Qx, [x^2 + 1])
    res = msolve(I)
    C, x = PolynomialRing(QQ, "x")
    @test res[2] == []
    @test res[1][1] == x^2+1
end
