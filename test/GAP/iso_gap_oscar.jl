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
   end
end
