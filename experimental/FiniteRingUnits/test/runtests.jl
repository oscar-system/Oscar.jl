function _test_map(f, D, C; morphism)
  @test domain(f) === D
  @test codomain(f) === C
  A = domain(f)
  for i in 1:10
    x = A isa FinGenAbGroup ? rand(A) : rand_pseudo(A; radius = 2)
    @test parent(f(x)) === C
  end

  if C isa FiniteRing
    a = rand(C)
    while !is_unit(a)
      a = rand(C)
    end
    b = rand(C)
    while !is_unit(b)
      b = rand(C)
    end
    if A isa FinGenAbGroup
      @test f\(a * b) == (f\a) + (f\b)
    end
  end

  if morphism
    if A isa FinGenAbGroup
      x, y = rand(A), rand(A)
      @test f(x + y) == f(x) * f(y)
    else
      x, y = rand_pseudo(A; radius = 2), rand_pseudo(A; radius = 2)
      @test f(x * y) == f(x) * f(y)
    end
  end
end

@testset "Unit groups" begin
  for (G, el) in [
                  #(symmetric_group(4), Dict(2 =>[2, 2, 4], 3 => [2, 2, 2, 6])),
                  (alternating_group(4), Dict(2 => [6], 3 => [6, 6])),
                  #(alternating_group(5), Dict(2 => [6], 3 => [2, 2, 24])),
                  (dihedral_group(8), Dict(2 => [2, 2, 4], 3 => [2, 2, 2, 2, 2])),
                  #(dihedral_group(16), Dict(2 => [2, 2, 2, 8], 3 => [2, 2, 2, 2, 2, 8])),
                  #(direct_product(cyclic_group(2), dihedral_group(8)), Dict(2 => [2, 2, 2, 2, 2, 2, 2, 4], 3 => [2, 2, 2, 2, 2, 2, 2, 2, 2, 2])),
                 ]
    for p in [2, 3]
      rngs = (GF(p)[G], finite_ring(GF(p)[G])[1])
      for R in rngs
        for (D, h) in [unit_group(R)]#, unit_group(DirectProductGroup, R)]
          _test_map(h, D, R; morphism = true)
        end
        Dab, hh = Oscar.FiniteRingUnits.abelianization_of_unit_group(R)
        @test elementary_divisors(Dab) == el[p]
        _test_map(hh, Dab, R; morphism = false)
      end
    end
  end
end

@testset "K1" begin
  for (G, el) in [#(symmetric_group(4), [2, 4]),
                  #(alternating_group(4), [6]),
                  #(alternating_group(5), [6]),
                  (dihedral_group(8), [2, 2, 4]),
                  #(dihedral_group(16), [2, 2, 2, 8]),
                  #(direct_product(cyclic_group(2), dihedral_group(8)), [2, 2, 2, 2, 2, 2, 2, 4]),
                 ]
    R, = finite_ring(GF(2)[G])
    D, h = Oscar.k1(R)
    @test elementary_divisors(D) == el
    _test_map(h, D, R; morphism = false)
  end

  let
    R, = finite_ring(matrix_algebra(GF(2), 2))
    D, f = Oscar.FiniteRingUnits.k1_semisimple_ring(R)
    @test order(D) == 1
    R, = finite_ring(direct_product(StructureConstantAlgebra(matrix_algebra(GF(2), 2))[1],
                                    StructureConstantAlgebra(matrix_algebra(GF(2), 2))[1])[1])
    D, f = Oscar.FiniteRingUnits.k1_semisimple_ring(R)
    _test_map(f, D, R; morphism = false)
    @test order(D) == 1
  end
end
