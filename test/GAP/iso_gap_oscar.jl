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
  # On the GAP side, consider
  # - fields created with `GF` (collections of `FFE`s):
  #   - of order up to 2^16 and larger fields,
  #   - with prescribed polynomial or with another polynomial.
  # - fields created with `AlgebraicExtension`.
  # On the Oscar side, consider fields of order less than 2^64 and larger ones.
#TODO: Eventually deal with GAP fields that are extensions of non-prime fields.
#      (See the lines below that are commented out.)
  @testset "with characteristic $p" for p in [ 5, 65537, next_prime(fmpz(2)^64) ]
    gap_p = GAP.Obj(p)
    F1 = GAP.Globals.GF(gap_p)
    R1 = GAP.Globals.PolynomialRing(F1);
#   F2 = GAP.Globals.GF(gap_p, 2)
#   R2 = GAP.Globals.PolynomialRing(F2);
    @testset "with degree $d" for d in 2:3
      pol1 = GAP.Globals.RandomPol(F1, d)
      while ! GAP.Globals.IsIrreducible(R1, pol1)
        pol1 = GAP.Globals.RandomPol(F1, d)
      end
#     pol2 = GAP.Globals.RandomPol(F2, d)
#     while ! GAP.Globals.IsIrreducible(R2, pol2)
#       pol2 = GAP.Globals.RandomPol(F2, d)
#     end
      @testset "with field $F" for F in [
             GAP.Globals.GF(gap_p, d),
             GAP.Globals.GF(F1, pol1),
#            GAP.Globals.AsField(F2, GAP.Globals.GF(gap_p, 2*d)),
             GAP.Globals.AlgebraicExtension(F1, pol1),
#            GAP.Globals.AlgebraicExtension(F2, pol2),
             ]
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
end

@testset "field of rationals" begin
  F = GAP.Globals.Rationals
  iso = Oscar.iso_gap_oscar(F)
  x = GAP.GapObj(2//3)
  y = 1
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
