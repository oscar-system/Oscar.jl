@testset "residue ring" begin
  @testset "order $n" for n in [ 2, 3, 65536, ZZRingElem(2)^64 ]
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
    @test_throws ErrorException preimage(iso, one(residue_ring(ZZ, n2)[1]))
  end
end

@testset "finite prime field" begin
  # On the GAP side, consider fields of order up to 2^16 and larger fields.
  # On the Oscar side, consider fields of order less than 2^64 and larger ones.
  @testset "with characteristic $p" for p in [ 5, 65537, next_prime(ZZRingElem(2)^64) ]
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
  @testset "with characteristic $p" for p in [5, 65537, next_prime(ZZRingElem(2)^64)]
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

  # Test the situation where a finite `AlgebraicExtension` field is chosen
  # on the GAP side.  (The test will fail if GAP will again run into the
  # problem of https://github.com/gap-system/gap/issues/4694.)
  R, x = polynomial_ring(GF(2))
  iso = Oscar.iso_oscar_gap(R)
  FG = GAP.Globals.AlgebraicExtension(
           GAP.Globals.LeftActingDomain(codomain(iso)), iso(x^2+x+1))
  iso = Oscar.iso_gap_oscar(FG)
  @test GAP.Globals.IsCanonicalBasisAlgebraicExtension(
            GAP.Globals.Basis(domain(iso)))
end

@testset "field of rationals, ring of integers" begin
  for (R, x, y) in [(GAP.Globals.Rationals, GapObj(2//3), 1),
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
    @test_throws ErrorException preimage(iso, cyclotomic_field(2)[2])
  end
end

@testset "number fields" begin  # includes quadratic fields
   @testset "subfield $F of a cyclotomic field" for F in
     [ GAP.evalstr( "Field( [ Sqrt(5) ] )" ),
       GAP.evalstr( "Field( [ Sqrt(-1) ] )" ),
       GAP.evalstr( "Field( [ Sqrt(6) ] )" ),
       GAP.evalstr( "Field( [ Sqrt(-6) ] )" ),
       GAP.evalstr( "Field( [ EC(19) ] )" ) ]  # not a quadratic field

     x = GAP.Globals.GeneratorsOfField(F)[1]
     y = GAP.Globals.One(F)
     iso = Oscar.iso_gap_oscar(F)
     @test iso === Oscar.iso_gap_oscar(F)  # test that everything gets cached
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

   R, x = polynomial_ring(QQ)
   iso = Oscar.iso_oscar_gap(R)
   QQG = GAP.Globals.LeftActingDomain(codomain(iso))

   @testset "`AlgebraicExtension` for polynomial $pol" for pol in
     [ x^2 - 5, x^2 + 3, x^3 - 2 ]

     F = GAP.Globals.AlgebraicExtension( QQG, iso(pol))
     x = GAP.Globals.GeneratorsOfField(F)[1]
     y = GAP.Globals.One(F)
     f = Oscar.iso_gap_oscar(F)
     @test f === Oscar.iso_gap_oscar(F)  # test that everything gets cached
     ox = f(x)
     oy = f(y)
     for i in 1:10
       xi = x^i
       oxi = f(xi)
       @test preimage(f, oxi) == xi
       @test oxi == ox^i
       @test oxi + oy == f(xi + y)
     end
   end

   @testset "`FieldByMatrices` isom. to field generated by $genstr" for genstr in
     ["1", "Sqrt(5)", "Sqrt(-7)"]
#    ["1", "Sqrt(5)", "Sqrt(-7)", "E(5)"]
#TODO: cyclotomic fields should support `GeneratorsOfAlgebraWithOne`

     F = GAP.Globals.Field(GAP.evalstr(genstr))
     matF = GAP.Globals.Image(GAP.Globals.IsomorphismMatrixAlgebra(F))
     mats = GAP.Globals.GeneratorsOfAlgebra(matF)
     F = GAP.Globals.FieldByMatrices(mats)
     x = GAP.Globals.GeneratorsOfField(F)[1]
     y = GAP.Globals.One(F)
#TODO: output of `FieldByMatrices` should support `IsFinite`
GAP.Globals.SetIsFinite(F, false)
     f = Oscar.iso_gap_oscar(F)
     @test f === Oscar.iso_gap_oscar(F)  # test that everything gets cached
     ox = f(x)
     oy = f(y)
     for i in 1:10
       xi = x^i
       oxi = f(xi)
       @test preimage(f, oxi) == xi
       @test oxi == ox^i
       @test oxi + oy == f(xi + y)
       @test preimage(f, oxi + oy) == xi + y
     end
   end
end

@testset "abelian closure" begin
  iso = Oscar.iso_gap_oscar(GAP.Globals.Cyclotomics)
  for N in [1, 2, 5, 15]
    x = GAP.Globals.E(N)
    y = iso(x)
    @test x == preimage(iso, y)
  end
  @test_throws ErrorException iso(GAP.Globals.Z(2))
  @test_throws ErrorException image(iso, GAP.Globals.Z(2))
  @test_throws ErrorException preimage(iso, cyclotomic_field(2)[2])
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
      @test_throws ErrorException preimage(iso, polynomial_ring(QQ, :y)[1]())
   end
end

@testset "multivariate polynomial rings" begin
   baserings = [GAP.Globals.Rationals,
                GAP.Globals.Integers,
                GAP.Globals.GF(2),
                GAP.Globals.GF(2, 3),
               ]
   @testset for R in baserings
      PR = GAP.Globals.PolynomialRing(R, 3)
      x, y, z = GAP.Globals.GeneratorsOfAlgebraWithOne(PR)
      iso = Oscar.iso_gap_oscar(PR)
      for pol in [zero(x), one(x), x, x^3+x+1, x*y+z+1]
         img = iso(pol)
         @test preimage(iso, img) == pol
      end
      @test_throws ErrorException iso(GAP.Globals.Z(2))
      @test_throws ErrorException image(iso, GAP.Globals.Z(2))
      @test_throws ErrorException preimage(iso, polynomial_ring(QQ, :y)[1]())
   end
end
