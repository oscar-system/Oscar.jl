@testset "computation and caching of homology" begin
  R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z])
  R3 = FreeMod(R, 3)
  v = sum(x*e for (x, e) in zip(gens(R), gens(R3)))
  K = Oscar.SimpleComplexWrapper(koszul_complex(v))

  KK = hom(K, K);

  tot = total_complex(KK)

  for i in upper_bound(tot):-1:lower_bound(tot)
    H, _ = homology(Oscar.underlying_complex(tot), 1, (i,))
    i > 0 && @test iszero(H)
    -3 <= i <= 0 && @test !iszero(H)
    -3 > i && @test iszero(H)
  end

  tot = simplify(tot)

  for i in upper_bound(tot):-1:lower_bound(tot)
    H, _ = homology(Oscar.underlying_complex(tot), 1, (i,))
    i > 0 && @test iszero(H)
    -3 <= i <= 0 && @test !iszero(H)
    -3 > i && @test iszero(H)
  end
end
