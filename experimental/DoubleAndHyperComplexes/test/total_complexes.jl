@testset "total complexes" begin
  R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])
  R1 = FreeMod(R, 1)
  Kx = Oscar.SimpleComplexWrapper(koszul_complex(x*R1[1]))
  Ky = Oscar.SimpleComplexWrapper(koszul_complex(y*R1[1]))
  Kz = Oscar.SimpleComplexWrapper(koszul_complex(z*R1[1]))
  Kw = Oscar.SimpleComplexWrapper(koszul_complex(w*R1[1]))

  K = hom(hom(Kx, Ky), hom(Kz, Kw))

  tot = total_complex(K);

  for i in upper_bound(tot):-1:lower_bound(tot)+2
    @test iszero(compose(map(tot, i), map(tot, i-1)))
  end
  @test !all(k->iszero(map(tot, k)), upper_bound(tot):-1:lower_bound(tot)+1)

  R4 = FreeMod(R, 4)
  v = sum(x*e for (x, e) in zip(gens(R), gens(R4)))
  KK = shift(Oscar.SimpleComplexWrapper(koszul_complex(v)), 2)

  phi = Oscar.lift_map(KK, tot, hom(KK[-2], tot[-2], [tot[-2][1]]), start_index=-2)

  @test all(k->is_isomorphism(phi[k]), -2:2)
end
