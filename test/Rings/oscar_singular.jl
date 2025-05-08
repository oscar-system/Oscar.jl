@testset "Oscar-Singular conversion" begin
  Qx, x = QQ[:x]
  K, a = number_field(x^3 + 2)
  S, = rational_function_field(K, "a")
  R1, = residue_ring(ZZ, 2)
  R2, = residue_ring(ZZ, ZZ(2)^100)
  FFrel = let # relative finite field
    _, x = GF(4)[:x]
    finite_field(x^3 + x + 1)[1]
  end
  Krel = let
    _, x = K[:x]
    number_field(x^2 - 3)[1]
  end
  FFt, = rational_function_field(GF(2), :t)
  FFtuv, = rational_function_field(GF(2), [:t, :u, :v])
  QQt, = rational_function_field(QQ, :t)
  QQtuv, = rational_function_field(QQ, [:t, :u, :v])

  test_rings = (ZZ, QQ, Qx,
                Nemo.Native.GF(2), GF(2), GF(2, 2), GF(next_prime(ZZ(2)^70)),
                GF(next_prime(ZZ(2)^70), 2), FFrel, abelian_closure(QQ)[1],
                K, Krel, FFt, QQt, FFtuv, QQtuv) 

  for R in test_rings
    f = Oscar.iso_oscar_singular_coeff_ring(R)
    @test domain(f) === R
    @test is_one(f(one(R)))
    @test is_zero(f(zero(R)))
    for i in 1:10
      a = R(rand(ZZ, -10:10))
      b = R(rand(ZZ, -10:10))
      @test f(a) + f(b) == f(a + b)
      @test f(a) * f(b) == f(a * b)
      @test preimage(f, f(a)) == a
    end

    Rx, (x, y) = polynomial_ring(R, [:x, :y])
    g = Oscar.iso_oscar_singular_poly_ring(Rx)
    @test domain(g) === Rx
    for i in 1:10
      a = R(rand(ZZ, -10:10))
      b = R(rand(ZZ, -10:10))
      e = rand(1:10)
      f = rand(1:10)
      h = a + x^e + b * y^f
      @test g(h)^2 == g(h^2)
      @test preimage(g, g(h)) == h
    end

    if R isa Field && !(R isa AbstractAlgebra.Generic.RationalFunctionField)
      Q, = quo(Rx, [x^2 - 1])
      x, y = gens(Q)
      g = Oscar.iso_oscar_singular_poly_ring(Q)
      @test domain(g) === Q
      for i in 1:10
        a = R(rand(ZZ, -10:10))
        b = R(rand(ZZ, -10:10))
        e = rand(1:10)
        f = rand(1:10)
        h = a + x^e + b * y^f
        @test g(h)^2 == g(h^2)
        @test preimage(g, g(h)) == h
      end
    end
  end

  let
    R, = rational_function_field(QQ, [:t, :u, :v])
    Rx, = R[:x, :y]
    f = Oscar.iso_oscar_singular_poly_ring(Rx)
    @test base_ring(codomain(f)) isa Singular.N_FField
  end
end
