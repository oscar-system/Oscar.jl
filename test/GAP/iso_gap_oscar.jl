@testset "residue ring" begin
  @testset "order $n" for n in [ 2, 3, 65536, fmpz(2)^64 ]
    gap_n = GAP.Obj(n)
    R = GAP.Globals.mod(GAP.Globals.Integers, gap_n)
    y = GAP.Globals.One(R)
    x = 2*y
    iso = Oscar.iso_gap_oscar(R)
    ox = iso(x)
    oy = iso(y)
    for i in 1:10
      xi = x^i
      oxi = iso(xi)
      @test preimage(iso, oxi) == xi
      @test oxi == ox^i
      @test oxi + oy == iso(xi + y)
    end
    n2 = n + 1

    x = GAP.Globals.ZmodnZObj(1, GAP.Obj(n2))
    @test_throws ErrorException iso(x)
    @test_throws ErrorException image(iso, x)
    @test_throws ErrorException preimage(iso, one(ResidueRing(ZZ, n2)))
  end
end

@testset "finite prime field" begin
  # On the GAP side, consider fields of order up to 2^16 and larger fields.
  # On the Oscar side, consider fields of order less than 2^64 and larger ones.
  @testset "with characteristic $p" for p in [ 5, 65537, next_prime(fmpz(2)^64) ]
    gap_p = GAP.Obj(p)
    F = GAP.Globals.GF(gap_p)
    x = GAP.Globals.Z(gap_p)
    y = GAP.Globals.One(F)
    iso = Oscar.iso_gap_oscar(F)
    ox = iso(x)
    oy = iso(y)
    for i in 1:10
      xi = x^i
      oxi = iso(xi)
      @test preimage(iso, oxi) == xi
      @test oxi == ox^i
      @test oxi + oy == iso(xi + y)
    end
    p2 = next_prime(p)
    @test_throws ErrorException iso(GAP.Globals.Z(GAP.Obj(p2)))
    @test_throws ErrorException image(iso, GAP.Globals.Z(GAP.Obj(p2)))
    @test_throws ErrorException preimage(iso, GF(p2)(1))
  end
end

@testset "finite non-prime field" begin
  # On the GAP side, consider only finite fields that consist of FFEs.
  # That is, avoid finite `AlgebraicExtension` fields and ignore
  # the defining polynomials of the GAP fields.
  # We support isomorphisms only from GAP fields that are extensions of
  # the prime field.
  # On the Oscar side, consider fields of order less than 2^64 and larger ones.
# @testset "with characteristic $p" for p in [5, 65537, next_prime(fmpz(2)^64)]
#T currently these the additional tests fail due to missing `embed_matrices`
  @testset "with characteristic $p" for p in [5]
    gap_p = GAP.Obj(p)
    @testset "with degree $d" for d in 2:3
      F = GAP.Globals.GF(gap_p, d)
      x = GAP.Globals.PrimitiveElement(F)
      y = GAP.Globals.One(F)
      iso = Oscar.iso_gap_oscar(F)
      ox = iso(x)
      oy = iso(y)
      for i in 1:10
        xi = x^i
        oxi = iso(xi)
        @test preimage(iso, oxi) == xi
        @test oxi == ox^i
        @test oxi + oy == iso(xi + y)
      end
      p2 = next_prime(p)
      o = GAP.Globals.One(GAP.Globals.GF(GAP.Obj(p2)))
      @test_throws ErrorException iso(o)
      @test_throws ErrorException image(iso, o)
      @test_throws ErrorException preimage(iso, GF(p2)(1))
    end
  end
end

@testset "field of rationals, ring of integers" begin
  for (R, x, y) in [(GAP.Globals.Rationals, GAP.GapObj(2//3), 1),
                    (GAP.Globals.Integers, 2, 3),
                   ]
    iso = Oscar.iso_gap_oscar(R)
    ox = iso(x)
    oy = iso(y)
    for i in 1:10
      xi = x^i
      oxi = iso(xi)
      @test preimage(iso, oxi) == xi
      @test oxi == ox^i
      @test oxi + oy == iso(xi + y)
    end
    @test_throws ErrorException preimage(iso, 1)
    @test_throws ErrorException iso(GAP.Globals.Z(2))
    @test_throws ErrorException image(iso, GAP.Globals.Z(2))
    @test_throws ErrorException preimage(iso, GF(2)(1))
  end
end

@testset "cyclotomic fields" begin
  @testset "with conductor $N" for N in [ 3, 4, 45 ]
    F = GAP.Globals.CF(N)
    x = GAP.Globals.E(N)
    y = GAP.Globals.One(F)
    iso = Oscar.iso_gap_oscar(F)
    ox = iso(x)
    oy = iso(y)
    for i in 1:10
      xi = x^i
      oxi = iso(xi)
      @test preimage(iso, oxi) == xi
      @test oxi == ox^i
      @test oxi + oy == iso(xi + y)
    end
    @test_throws ErrorException iso(GAP.Globals.Z(2))
    @test_throws ErrorException image(iso, GAP.Globals.Z(2))
    @test_throws ErrorException preimage(iso, CyclotomicField(2)[2])
  end
end

@testset "univariate polynomial rings" begin
   baserings = [GAP.Globals.Rationals,
                GAP.Globals.Integers,
                GAP.Globals.GF(2),
                GAP.Globals.GF(2, 3),
               ]
   @testset for R in baserings
      PR = GAP.Globals.PolynomialRing(R)
      x = GAP.Globals.Indeterminate(R)
      iso = Oscar.iso_gap_oscar(PR)
      for pol in [zero(x), one(x), x, x^3+x+1]
         img = iso(pol)
         @test preimage(iso, img) == pol
      end
      @test_throws ErrorException iso(GAP.Globals.Z(2))
      @test_throws ErrorException image(iso, GAP.Globals.Z(2))
      @test_throws ErrorException preimage(iso, PolynomialRing(QQ, "y")[1]())
   end
end
