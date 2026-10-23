@testset "FunctionField" begin
  for k in [QQ, GF(3)]
    tests = []

    # rational function field
    F, a = rational_function_field(k, "a")
    push!(tests, (F, a, ["a"]))

    # univariate fraction_field
    F = fraction_field(k[:a][1])
    a = F(gen(base_ring(F)))
    push!(tests, (F, a, ["a"]))

    # multivariate fraction_field
    F = fraction_field(k[:a1, :a2][1])
    a = F(base_ring(F)[1])
    push!(tests, (F, a, ["a1", "a2"]))

    for (F, a, symb) in tests
      coeff_iso = Oscar.iso_oscar_singular_coeff_ring(F)
      @test codomain(coeff_iso) isa Singular.N_FField
      @test preimage(coeff_iso, coeff_iso(a)) == a
      @test preimage(coeff_iso, coeff_iso(inv(a))) == inv(a)

      R, (x, y) = polynomial_ring(F, [:x, :y])
      iso = Oscar.iso_oscar_singular_poly_ring(R)
      @test codomain(iso) isa Singular.PolyRing{<:Singular.n_transExt}

      @test k(F(1)) isa elem_type(k)
      @test_throws InexactError k(a)
    end
  end
end
