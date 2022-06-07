@testset "residue rings" begin
   @testset for n in [2, 3, 6]
      for R in [ResidueRing(ZZ, n), ResidueRing(ZZ, fmpz(n))]
         f = Oscar.iso_oscar_gap(R)
         a = one(R)
         b = -one(R)
         @test f(a*b) == f(a)*f(b)
         @test f(a-b) == f(a)-f(b)
         C = codomain(f)
         for a in C
            for b in C
               @test preimage(f, a*b) == preimage(f, a)*preimage(f, b)
               @test preimage(f, a-b) == preimage(f, a)-preimage(f, b)
            end
         end
         n2 = n + 1
         one2 = one(ResidueRing(ZZ, n2))
         @test_throws ErrorException f(one2)
         @test_throws ErrorException image(f, one2)
         @test_throws ErrorException preimage(f, GAP.Globals.ZmodnZObj(1, GAP.Obj(n2)))
      end
   end

   n = fmpz(2)^100
   R = ResidueRing(ZZ, n)
   f = Oscar.iso_oscar_gap(R)
   a = -one(R)
   @test f(a) == -f(one(R))
   C = codomain(f)
   a = -GAP.Globals.One(C)
   @test preimage(f, a) == -preimage(f, GAP.Globals.One(C))
end

@testset "finite fields" begin
   @testset for p in [2, 3]
      for F in [Nemo.GaloisField(UInt(p)), Nemo.GaloisFmpzField(fmpz(p))]
         f = Oscar.iso_oscar_gap(F)
         for a in F
            for b in F
               @test f(a*b) == f(a)*f(b)
               @test f(a-b) == f(a)-f(b)
            end
         end
         C = codomain(f)
         for a in C
            for b in C
               @test preimage(f, a*b) == preimage(f, a)*preimage(f, b)
               @test preimage(f, a-b) == preimage(f, a)-preimage(f, b)
            end
         end
         p2 = next_prime(p)
         @test_throws ErrorException f(GF(p2)(1))
         @test_throws ErrorException image(f, GF(p2)(1))
         @test_throws ErrorException preimage(f, GAP.Globals.Z(GAP.Obj(p2)))
      end
   end

   @testset for (p,d) in [(2, 1), (5, 1), (2, 4), (3, 3)]
      for F in [FqFiniteField(fmpz(p),d,:z), FqDefaultFiniteField(fmpz(p),d,:z)]
         f = Oscar.iso_oscar_gap(F)
         g = elm -> map_entries(f, elm)
         for a in F
            for b in F
               @test f(a*b) == f(a)*f(b)
               @test f(a-b) == f(a)-f(b)
            end
         end
         C = codomain(f)
         for a in C
            for b in C
               @test preimage(f, a*b) == preimage(f, a)*preimage(f, b)
               @test preimage(f, a-b) == preimage(f, a)-preimage(f, b)
            end
         end
         G = GL(4,F)
         for a in gens(G)
            for b in gens(G)
               @test g(a.elm*b.elm) == g(a.elm)*g(b.elm)
               @test g(a.elm-b.elm) == g(a.elm)-g(b.elm)
            end
         end
         p2 = next_prime(p)
         @test_throws ErrorException f(GF(p2)(1))
         @test_throws ErrorException image(f, GF(p2)(1))
         @test_throws ErrorException preimage(f, GAP.Globals.Z(GAP.Obj(p2)))
      end
   end
end

@testset "a large non-prime field (FqNmodFiniteField)" begin
   # The defining polynomial of the Oscar field is not a Conway polynomial,
   # the polynomial of the GAP field is a Conway polynomial,
   # thus we need an intermediate field on the Oscar side.
   p = next_prime(10^6)
   F = GF(p, 2)
   f = Oscar.iso_oscar_gap(F)
   for x in [ F(3), gen(F) ]
      a = f(x)
      @test preimage(f, a) == x
   end
   @test GAP.Globals.DefiningPolynomial(codomain(f)) ==
         GAP.Globals.ConwayPolynomial(p, 2)
   @test F.is_conway == 0
   p2 = next_prime(p)
   @test_throws ErrorException f(GF(p2)(1))
   @test_throws ErrorException image(f, GF(p2)(1))
   @test_throws ErrorException preimage(f, GAP.Globals.Z(GAP.Obj(p2)))
end

@testset "another large non-prime field (FqNmodFiniteField)" begin
   # GAP's `GF(p, d)` throws an error if the Conway polynomial in question
   # is neither known nor cheap to compute.
   # Here we can translate from Oscar to GAP by choosing the
   # defining polynomial of the Oscar field for constructing the GAP field.
   p = 1031
   d = 10
   F = GF(p, d)

   # Check whether the current example is still relevant.
   @test ! GAP.Globals.IsCheapConwayPolynomial(p, d)

   f = Oscar.iso_oscar_gap(F)
   for x in [ F(3), gen(F) ]
      a = f(x)
      @test preimage(f, a) == x
   end
   p2 = next_prime(p)
   @test_throws ErrorException f(GF(p2)(1))
   @test_throws ErrorException image(f, GF(p2)(1))
   @test_throws ErrorException preimage(f, GAP.Globals.Z(GAP.Obj(p2)))
end

@testset "field of rationals, ring of integers" begin
  for (R, x, y) in [(QQ, QQ(2//3), QQ(1)), (ZZ, ZZ(2), ZZ(1))]
    iso = Oscar.iso_oscar_gap(R)
    ox = iso(x)
    oy = iso(y)
    for i in 1:10
      xi = x^i
      oxi = iso(xi)
      @test preimage(iso, oxi) == xi
      @test oxi == ox^i
      @test oxi + oy == iso(xi + y)
    end
    @test_throws ErrorException iso(1)
    @test_throws ErrorException image(iso, 1)
    @test_throws ErrorException iso(GF(2)(1))
    @test_throws ErrorException image(iso, GF(2)(1))
    @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
  end
end

@testset "cyclotomic fields" begin
   # for computing random elements of the fields in question
   my_rand_bits(F::FlintRationalField, b::Int) = rand_bits(F, b)
   my_rand_bits(F::AnticNumberField, b::Int) = F([rand_bits(QQ, b) for i in 1:degree(F)])

   fields = Any[CyclotomicField(n) for n in [1, 3, 4, 5, 8, 15, 45]]
   push!(fields, (QQ, 1))

   @testset for (F, z) in fields
      f = Oscar.iso_oscar_gap(F)
      g = elm -> map_entries(f, elm)
      for i in 1:10
         a = my_rand_bits(F, 5)
         for j in 1:10
            b = my_rand_bits(F, 5)
            @test f(a*b) == f(a)*f(b)
            @test f(a - b) == f(a) - f(b)
         end
      end
      @test_throws ErrorException f(CyclotomicField(2)[2])
      @test_throws ErrorException image(f, CyclotomicField(2)[2])
      @test_throws ErrorException preimage(f, GAP.Globals.Z(2))
   end

   K, a = CyclotomicField(10, "a")
   phi = Oscar.iso_oscar_gap(K)
   @test phi(-a^3 + a^2 + 1) == GAP.evalstr("-E(5)^2-E(5)^3")
end

@testset "abelian closure" begin
  F, z = abelian_closure(QQ)
  iso = Oscar.iso_oscar_gap(F)
  for N in [1, 2, 5, 15]
    x = z(N)
    y = iso(x)
    @test x == preimage(iso, y)
  end
  @test_throws ErrorException iso(CyclotomicField(2)[2])
  @test_throws ErrorException image(iso, CyclotomicField(2)[2])
  @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
end

@testset "univariate polynomial rings" begin
   baserings = [QQ,                           # yields `FmpqPolyRing`
                ZZ,                           # yields `FmpzPolyRing`
                GF(2,2),                      # yields `FqNmodPolyRing`
#               FqDefaultFiniteField(fmpz(2),3,:x), # yields `FqDefaultPolyRing`
#TODO: This case fails due to a problem in AA/Nemo,
#      see https://github.com/Nemocas/Nemo.jl/issues/1307
                FqFiniteField(fmpz(2),2,:z),  # yields `FqPolyRing`
                GF(fmpz(2)),                  # yields `GFPFmpzPolyRing`
                GF(2),                        # yields `GFPPolyRing`
                Nemo.NmodRing(UInt64(6)),     # yields `NmodPolyRing`
                Nemo.FmpzModRing(fmpz(6)),    # yields `FmpzModPolyRing`
               ]
#TODO: How to get `AbstractAlgebra.Generic.PolyRing`?
   @testset for R in baserings
      PR, x = PolynomialRing(R, "x")
      iso = Oscar.iso_oscar_gap(PR)
      for pol in [zero(x), one(x), x, x^3+x+1]
         img = iso(pol)
         @test preimage(iso, img) == pol
      end
      m = matrix([x x; x x])
      @test map_entries(inv(iso), map_entries(iso, m)) == m
      @test_throws ErrorException iso(PolynomialRing(R, "y")[1]())
      @test_throws ErrorException image(iso, PolynomialRing(R, "y")[1]())
      @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
   end
end

@testset "multivariate polynomial rings" begin
   baserings = [QQ,                           # yields `FmpqMPolyRing`
                ZZ,                           # yields `FmpzMPolyRing`
                GF(2,2),                      # yields `FqNmodMPolyRing`
                GF(fmpz(2)),                  # yields `AbstractAlgebra.Generic.MPolyRing{gfp_fmpz_elem}`
                GF(2),                        # yields `GFPMPolyRing`
                Nemo.NmodRing(UInt64(6)),     # yields `NmodMPolyRing`
               ]
#TODO: How to get `GFPFmpzMPolyRing`, `FqDefaultMPolyRing`?
   @testset for R in baserings
      PR, (x,y,z) = PolynomialRing(R, 3)
      iso = Oscar.iso_oscar_gap(PR)
      for pol in [zero(x), one(x), x, x^2*y + y*z^3 + x*y*z + 1]
         img = iso(pol)
         @test preimage(iso, img) == pol
      end
      m = matrix([x x; y z])
      @test map_entries(inv(iso), map_entries(iso, m)) == m
      @test_throws ErrorException iso(PolynomialRing(R, ["y"])[1]())
      @test_throws ErrorException image(iso, PolynomialRing(R, ["y"])[1]())
      @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
   end
end
