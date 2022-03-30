@testset "FunctionField" begin
  for (k, kk) in [(QQ, Singular.QQ), (GF(3), Singular.Fp(3))]
    tests = []
    
    # RationalFunctionField
    F, a = RationalFunctionField(k, "a")
    push!(tests, (F, a, ["a"]))
    
    # univariate FractionField
    F = FractionField(k["a"][1])
    a = F(gen(base_ring(F)))
    push!(tests, (F, a, ["a"]))
    
    # multivariate FractionField
    F = FractionField(k["a1", "a2"][1])
    a = F(base_ring(F)[1])
    push!(tests, (F, a, ["a1", "a2"]))
    
    for (F, a, symb) in tests
      F1, (a1,) = Singular.FunctionField(kk, symb)
      F2, (a2,) = Singular.FunctionField(kk, vcat(symb, ["b"])) # one more gen
      F3, (a3,) = Singular.FunctionField(Singular.Fp(5), symb) # different char
      
      @test singular_coeff_ring(F) === F1
      @test F(a1) == a
      @test F1(a) == a1
      
      @test_throws ArgumentError F(a2) # wrong ngens
      @test_throws ArgumentError F2(a) # wrong ngens
      
      @test_throws ArgumentError F(a3) # wrong char
      @test_throws ArgumentError F3(a) # wrong char
      
      R, (x, y) = PolynomialRing(F, ["x", "y"])
      @test singular_poly_ring(R) isa Singular.PolyRing{Singular.n_transExt}

      @test k(F(1)) isa elem_type(k)
      @test_throws InexactError k(a)
    end
  end
end
