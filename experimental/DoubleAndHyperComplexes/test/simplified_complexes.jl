@testset "simplified complexes" begin
  # We create two Koszul complexes, one which can not be simplified
  # and another one which is homotopy equivalent to the zero complex.
  n = 4
  d1 = 3
  R, x = polynomial_ring(QQ, n)
  Rn = FreeMod(R, n)
  v = sum(x^d1*e for (x, e) in zip(gens(R), gens(Rn)))
  K = Oscar.SimpleComplexWrapper(koszul_complex(v))

  KK = simplify(K) # This should change nothing, because there are 
                   # no units in the matrices.
  @test rank(KK[0]) == 1
  @test rank(KK[2]) == 6
  @test rank(KK[1]) == 4
  @test matrix(map(K, 1, (1,))) == matrix(map(KK, 1, (1,)))
  @test matrix(map(K, 1, (2,))) == matrix(map(KK, 1, (2,)))
  @test matrix(map(K, 1, (3,))) == matrix(map(KK, 1, (3,)))

  F = FreeMod(R, n+1)
  v = sum(x^d1*e for (x, e) in zip(gens(R), gens(F)[1:n]))
  v = v + F[n+1] # We add one more entry with a unit.

  K = Oscar.SimpleComplexWrapper(koszul_complex(v))

  KK = simplify(K) # This should contract everything to the zero complex.
  all(k->iszero(KK[k]), -1:n+1)
  @test !can_compute_index(KK, -2)
end

@testset "determinantal surface" begin
  # We use the Hilbert-Burch theorem here. 
  # The matrix A below appears again in the resolution 
  # of the quotient. We localize at the entry x-1 of 
  # the matrix and simplification should turn the 
  # complex into a Koszul complex of length 2.
  R, (x, y, z, w) = polynomial_ring(QQ, [:x, :y, :z, :w])

  A = R[x y z; y z x-1]

  I = ideal(R, minors(A, 2))

  Q, _ = quo(R, I)

  L, _ = localization(R, powers_of_element(x-1))

  L1 = FreeMod(L, 1)
  IL1, inc = L(I)*L1

  M = cokernel(inc)

  c, _ = free_resolution(Oscar.SimpleFreeResolution, M)

  cc = simplify(c) # Should be a Koszul complex of length 2.

  @test rank(c[0]) == 1
  @test rank(cc[0]) == 1
  @test rank(c[1]) == 3
  @test rank(cc[1]) == 2
  @test rank(c[2]) == 2
  @test rank(cc[2]) == 1

  for i in 2:-1:1
    orig = map(c, i)
    new_map = map(cc, i)
    # Check that all the squares commute.
    # We can get the maps relating cc with c from cc 
    # using map_from/to_original_complex.
    @test compose(orig, map_from_original_complex(cc)[i-1]) == compose(map_from_original_complex(cc)[i], new_map)
    @test compose(new_map, map_to_original_complex(cc)[i-1]) == compose(map_to_original_complex(cc)[i], orig)
  end
end

@testset "homotopies for simplified complexes" begin
  R, (x, y, z, w) = graded_polynomial_ring(QQ, [:x, :y, :z, :w])
  a = [x, one(R), y, z, -2*one(R), w]
  k = Oscar.HomogKoszulComplex(R, a)
  s = simplify(k);
  @test_throws ErrorException Oscar.homotopy_map(s, 3)
  s = simplify(k; with_homotopy_maps=true);

  for k in [rand(0:length(a)) for _ in 1:3*length(a)]
    s[k]
  end


  from = Oscar.map_from_original_complex(s);
  to = Oscar.map_to_original_complex(s);

  h0 = Oscar.homotopy_map(s, 0)
  @test compose(h0, map(k, 1)) == id_hom(k[0]) - compose(from[0], to[0])

  for i in 1:length(a)-1
    s[i]
    h0 = Oscar.homotopy_map(s, i)
    h1 = Oscar.homotopy_map(s, i-1)
    @test compose(h0, map(k, i+1)) + compose(map(k, i), h1) == id_hom(k[i]) - compose(from[i], to[i])
  end

  h1 = Oscar.homotopy_map(s, length(a)-1)
  @test compose(map(k, length(a)), h1) == id_hom(k[length(a)]) - compose(from[length(a)], to[length(a)])
end
