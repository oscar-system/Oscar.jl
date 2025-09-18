@testset "residue rings" begin
   @testset for n in [2, 3, 6]
      for R in [residue_ring(ZZ, n)[1], residue_ring(ZZ, ZZRingElem(n))[1]]
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
         one2 = one(residue_ring(ZZ, n2)[1])
         @test_throws ErrorException f(one2)
         @test_throws ErrorException image(f, one2)
         @test_throws ErrorException preimage(f, GAP.Globals.ZmodnZObj(1, GAP.Obj(n2)))
      end
   end

   n = ZZRingElem(2)^100
   R = residue_ring(ZZ, n)[1]
   f = Oscar.iso_oscar_gap(R)
   a = -one(R)
   @test f(a) == -f(one(R))
   C = codomain(f)
   a = -GAP.Globals.One(C)
   @test preimage(f, a) == -preimage(f, GAP.Globals.One(C))
end

@testset "finite fields" begin
   @testset for p in [2, 3]
      for F in [fpField(UInt(p)), FpField(ZZRingElem(p))]
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

   p = 257  # GAP regards the Conway polynomial for `GF(257, 1)` as not cheap.
   @testset for F in [fpField(UInt(257)), FpField(ZZRingElem(257))]
      f = Oscar.iso_oscar_gap(F)
      oO = one(F)
      oG = f(oO)
      for a in F
         @test f(a*a) == f(a)*f(a)
         @test f(a-oO) == f(a)-oG
      end
      C = codomain(f)
      for a in C
         @test preimage(f, a*a) == preimage(f, a)*preimage(f, a)
      end
   end

   @testset for (p,d) in [(2, 1), (5, 1), (2, 4), (3, 3)]
      for F in [FqPolyRepField(ZZRingElem(p),d,:z), FqField(ZZRingElem(p),d,:z)]
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
               @test g(matrix(a)*matrix(b)) == g(matrix(a))*g(matrix(b))
               @test g(matrix(a)-matrix(b)) == g(matrix(a))-g(matrix(b))
            end
         end
         p2 = next_prime(p)
         @test_throws ErrorException f(GF(p2)(1))
         @test_throws ErrorException image(f, GF(p2)(1))
         @test_throws ErrorException preimage(f, GAP.Globals.Z(GAP.Obj(p2)))
      end
   end
end

@testset "a large non-prime field (fqPolyRepField)" begin
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
   p2 = next_prime(p)
   @test_throws ErrorException f(GF(p2)(1))
   @test_throws ErrorException image(f, GF(p2)(1))
   @test_throws ErrorException preimage(f, GAP.Globals.Z(GAP.Obj(p2)))
end

@testset "another large non-prime field (fqPolyRepField)" begin
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
   my_rand_bits(F::QQField, b::Int) = rand_bits(F, b)
   my_rand_bits(F::AbsSimpleNumField, b::Int) = F([rand_bits(QQ, b) for i in 1:degree(F)])

   fields = Any[cyclotomic_field(n) for n in [1, 3, 4, 5, 8, 15, 45]]
   push!(fields, (QQ, 1))

   @testset for (F, z) in fields
      f = Oscar.iso_oscar_gap(F)
#     g = elm -> map_entries(f, elm)
#TODO: add some tests for mapping matrices over F
      for i in 1:10
         a = my_rand_bits(F, 5)
         for j in 1:10
            b = my_rand_bits(F, 5)
            @test f(a*b) == f(a)*f(b)
            @test f(a - b) == f(a) - f(b)
            @test preimage(f, f(a)) == a
         end
      end
      @test_throws ErrorException f(cyclotomic_field(2)[2])
      @test_throws ErrorException image(f, cyclotomic_field(2)[2])
      @test_throws ErrorException preimage(f, GAP.Globals.Z(2))
   end

   K, a = cyclotomic_field(10, "a")
   phi = Oscar.iso_oscar_gap(K)
   @test phi(-a^3 + a^2 + 1) == GAP.evalstr("-E(5)^2-E(5)^3")
end

@testset "quadratic number fields" begin
   # for computing random elements of the fields in question
   my_rand_bits(F::QQField, b::Int) = rand_bits(F, b)
   my_rand_bits(F::AbsSimpleNumField, b::Int) = F([rand_bits(QQ, b) for i in 1:degree(F)])

   @testset for N in [ 5, -3, 12, -8 ]
      F, z = quadratic_field(N)
      f = Oscar.iso_oscar_gap(F)
      @test f === Oscar.iso_oscar_gap(F)  # test that everything gets cached
      for i in 1:10
         a = my_rand_bits(F, 5)
         for j in 1:10
            b = my_rand_bits(F, 5)
            @test f(a*b) == f(a)*f(b)
            @test f(a - b) == f(a) - f(b)
            @test preimage(f, f(a)) == a
         end
      end

      m = matrix([z z; z z])
      @test map_entries(inv(f), map_entries(f, m)) == m
   end
end

@testset "number fields" begin
   # for computing random elements of the fields in question
   my_rand_bits(F::QQField, b::Int) = rand_bits(F, b)
   my_rand_bits(F::NumField, b::Int) = F([my_rand_bits(base_field(F), b) for i in 1:degree(F)])

   # absolute number fields
   R, x = polynomial_ring(QQ, :x)
   pols = [ x - 1, x^2 - 5, x^2 + 3, x^3 - 2,  # simple
            [x^2 - 2, x^2 + 1] ]        # non-simple
   fields = Any[number_field(pol)[1] for pol in pols]

   # non-absolute number fields
   F1, _ = number_field(x^2-2)
   R1, x1 = polynomial_ring(F1, :x)
   F2, _ = number_field(x1^2-3)
   push!(fields, F2)                                 # simple
   push!(fields, number_field([x1^2-3, x1^2+1])[1])  # non-simple

   # and a simple three-step construction (in order to exercise recursion)
   R2, x2 = polynomial_ring(F2, :x)
   push!(fields, number_field(x2^2+1)[1])

   @testset for F in fields
      f = Oscar.iso_oscar_gap(F)
      @test f === Oscar.iso_oscar_gap(F)  # test that everything gets cached
      for i in 1:10
         a = my_rand_bits(F, 5)
         for j in 1:10
            b = my_rand_bits(F, 5)
            @test f(a*b) == f(a)*f(b)
            @test f(a - b) == f(a) - f(b)
            @test preimage(f, f(a)) == a
         end
      end
   end

   # an application
   K = fields[5]
   a, b = gens(K)
   M1 = 1/a*matrix(K, [1 1; 1 -1])
   M2 = matrix(K, [1 0 ; 0 b])
   G = matrix_group(M1, M2)
   @test small_group_identification(G) == (192, 963)
end

@testset "abelian closure" begin
  F, z = abelian_closure(QQ)
  iso = Oscar.iso_oscar_gap(F)
  for N in [1, 2, 5, 15]
    x = z(N)
    y = iso(x)
    @test x == preimage(iso, y)
  end
  @test_throws ErrorException iso(cyclotomic_field(2)[2])
  @test_throws ErrorException image(iso, cyclotomic_field(2)[2])
  @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
end

@testset "univariate polynomial rings" begin
   baserings = [QQ,                           # yields `QQPolyRing`
                ZZ,                           # yields `ZZPolyRing`
                GF(2,3),                      # yields `FqPolyRing`
                Nemo.Native.GF(2,2),          # yields `fqPolyRepPolyRing`
                FqPolyRepField(ZZRingElem(2),2,:z),  # yields `FqPolyRepPolyRing`
                Nemo.Native.GF(ZZRingElem(2)),# yields `FpPolyRing`
                Nemo.Native.GF(2),            # yields `fpPolyRing`
                zzModRing(UInt64(6)),         # yields `zzModPolyRing`
                ZZModRing(ZZRingElem(6)),     # yields `ZZModPolyRing`
               ]
#TODO: How to get `AbstractAlgebra.Generic.PolyRing`?
   @testset for R in baserings
      PR, x = polynomial_ring(R, :x)
      iso = Oscar.iso_oscar_gap(PR)
      for pol in [zero(x), one(x), x, x^3+x+1]
         img = iso(pol)
         @test preimage(iso, img) == pol
      end
      m = matrix([x x; x x])
      @test map_entries(inv(iso), map_entries(iso, m)) == m
      @test_throws ErrorException iso(polynomial_ring(R, :y)[1]())
      @test_throws ErrorException image(iso, polynomial_ring(R, :y)[1]())
      @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
   end
end

@testset "multivariate polynomial rings" begin
   baserings = [QQ,                           # yields `QQMPolyRing`
                ZZ,                           # yields `ZZMPolyRing`
                GF(2,3),                      # yields `FqMPolyRing`
                Nemo.Native.GF(2,2),          # yields `fqPolyRepMPolyRing`
                FqPolyRepField(ZZRingElem(2),2,:z),  # yields `AbstractAlgebra.Generic.MPolyRing{FqPolyRepFieldElem}`
                Nemo.Native.GF(ZZRingElem(2)),# yields `FpMPolyRing`
                Nemo.Native.GF(2),            # yields `fpMPolyRing`
                zzModRing(UInt64(6)),         # yields `zzModMPolyRing`
               ]
#TODO: How to get `FpMPolyRing`, `FqMPolyRing`?
   @testset for R in baserings
      PR, (x,y,z) = polynomial_ring(R, 3)
      iso = Oscar.iso_oscar_gap(PR)
      for pol in [zero(x), one(x), x, x^2*y + y*z^3 + x*y*z + 1]
         img = iso(pol)
         @test preimage(iso, img) == pol
      end
      m = matrix([x x; y z])
      @test map_entries(inv(iso), map_entries(iso, m)) == m
      @test_throws ErrorException iso(polynomial_ring(R, [:y])[1]())
      @test_throws ErrorException image(iso, polynomial_ring(R, [:y])[1]())
      @test_throws ErrorException preimage(iso, GAP.Globals.Z(2))
   end
end
