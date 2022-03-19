@testset "MPolyQuo.graded" begin
  R, (x, y) = grade(PolynomialRing(QQ, ["x", "y"])[1]);
  Q = quo(R, ideal([x^2, y]))[1];
  h = homogeneous_components(Q[1])
  @test valtype(h) === elem_type(Q)
end

function _random_poly(RR, n)
  pols = elem_type(RR)[]
  for t in 1:n
    f = MPolyBuildCtx(RR)
    g = MPolyBuildCtx(RR.R)
    for z in 1:3
      e = rand(0:6, ngens(RR.R))
      r = base_ring(RR.R)(rand(0:22))
      push_term!(f, r, e)
      push_term!(g, r, e)
    end
    f = finish(f)
    @test f == RR(finish(g))
    push!(pols, RR(finish(g)))
  end
  return pols
end

function _homogeneous_polys(polys::Vector{<:MPolyElem})
  R = parent(polys[1])
  D = []
  _monomials = [zero(R)]
  for i = 1:4
    push!(D, homogeneous_components(polys[i]))
    _monomials = append!(_monomials, monomials(polys[i]))
  end
  hom_polys = []
  for deg in unique([iszero(mon) ? id(grading_group(R)) : degree(mon) for mon = _monomials])
    g = zero(R)
    for i = 1:4
      if haskey(D[i], deg)
        temp = get(D[i], deg, 'x')
        @test homogeneous_component(polys[i], deg) == temp
        g += temp
      end
    end
    @test ishomogeneous(g)
    push!(hom_polys, g)
  end
  return hom_polys
end

@testset "mpoly-graded" begin

    Qx, (x,y,z) = PolynomialRing(QQ, ["x", "y", "z"])
    t = gen(Hecke.Globals.Qx)
    k1 , l= number_field(t^5+t^2+2)
    NFx = PolynomialRing(k1, ["x", "y", "z"])[1]
    k2 = Nemo.GF(23)
    GFx = PolynomialRing(k2, ["x", "y", "z"])[1]
    # TODO explain why test fails if Nemo.ResidueRing(ZZ,17) => Nemo.GF(17)
    RNmodx=PolynomialRing(Nemo.ResidueRing(ZZ,17), :x => 1:2)[1]
    Rings= [Qx, NFx, GFx, RNmodx]

    A = abelian_group([4 3 0 1;
                       0 0 0 3])

    GrpElems = [A(fmpz[1,0,2,1]), A(fmpz[1,1,0,0]), A(fmpz[2,2,1,0])]

    T, _ = grade(Qx)
    @test sprint(show, "text/plain", T(x)//T(1)) isa String

    v = [1,1,2]

    for R in Rings
      decorated_rings = [decorate(R)[1],
                         grade(R, [v[i] for i=1:ngens(R)])[1],
                         filtrate(R, [v[i] for i=1:ngens(R)])[1],
                         filtrate(R, [GrpElems[i] for i =1:ngens(R)], (x, y) -> x[1]+x[2]+x[3]+x[4] < y[1]+y[2]+y[3]+y[4])[1],
                         grade(R, [GrpElems[i] for i =1:ngens(R)])[1]]
    end

    for R in Rings
      decorated_rings = [decorate(R)[1],
                         grade(R, [v[i] for i=1:ngens(R)])[1],
                         filtrate(R, [v[i] for i=1:ngens(R)])[1],
                         filtrate(R, [GrpElems[i] for i =1:ngens(R)], (x, y) -> x[1]+x[2]+x[3]+x[4] < y[1]+y[2]+y[3]+y[4])[1],
                         grade(R, [GrpElems[i] for i =1:ngens(R)])[1]]

      graded_rings = decorated_rings[[2, 5]]
      filtered_rings =  decorated_rings[[1, 3, 4]]

      d_Elems = IdDict(decorated_rings[1] => 4*decorated_rings[1].D[1],
                       decorated_rings[2] => 5*decorated_rings[2].D[1],
                       decorated_rings[3] => 3*decorated_rings[3].D[1],
                       decorated_rings[4] => decorated_rings[4].D([2,2,1,0]),
                       decorated_rings[5] => decorated_rings[5].D([2,2,1,0]))

      dims = IdDict(zip(decorated_rings, [15, 12, 6, 1, 1]))

      # In RR, we take the quotient by the homogeneous component corresponding to d_Elems[RR]
      # The dimension should be dims[RR]

      for RR in decorated_rings
        @test (RR in graded_rings) == Oscar.isgraded(RR)
        @test (RR in filtered_rings) == Oscar.isfiltered(RR)
        polys = _random_poly(RR, 4) # create 4 random polynomials
        @test ngens(RR) == length(gens(RR))
        @test gen(RR, 1) == RR[1]

        @test one(RR) == RR(1)
        @test isone(one(RR))
        @test zero(RR) == RR(0)
        @test iszero(zero(RR))
        @test !iszero(one(RR))
        @test !isone(zero(RR))
        @test divexact(one(RR), one(RR)) == one(RR)

        @test (polys[1] + polys[2])^2 == polys[1]^2 + 2*polys[1]*polys[2] + polys[2]^2
        @test (polys[3] - polys[4])^2 == polys[3]^2 + 2*(-polys[3])*polys[4] + polys[4]^2
        @test polys[2] * (polys[3] + polys[4]) == Oscar.addeq!(Oscar.mul!(polys[1], polys[2], polys[3]), Oscar.mul!(polys[1], polys[2], polys[4]))
        @test parent(polys[1]) === RR

        for k in 1:length(polys[4])
          @test coeff(polys[4],k) * Oscar.monomial(polys[4], k) ==
                finish(push_term!(MPolyBuildCtx(RR), collect(Oscar.MPolyCoeffs(polys[4]))[k], collect(Oscar.MPolyExponentVectors(polys[4]))[k]))
        end

        hom_polys = _homogeneous_polys(polys)
        I = ideal(hom_polys)
        R_quo = Oscar.MPolyQuo(RR, I)
        @test base_ring(R_quo) == RR
        @test modulus(R_quo) == I
        f = R_quo(polys[2])
        D = homogeneous_components(f)
        for deg in [degree(R_quo(mon)) for mon  = Oscar.monomials(f.f)]
          h = get(D, deg, 'x')
          @test ishomogeneous(R_quo(h))
          @test h == R_quo(homogeneous_component(f, deg))
        end

        @test Oscar.isfiltered(R_quo) == Oscar.isfiltered(RR)
        @test Oscar.isgraded(R_quo) == Oscar.isgraded(RR)

        @test grading_group(R_quo) == grading_group(RR)

        d_Elem = d_Elems[RR]

        dim_test = dims[RR]

        if base_ring(R) isa AbstractAlgebra.Field
	  if isfree(grading_group(RR))
             grp_elem = d_Elem
             H = homogeneous_component(RR, grp_elem)
             @test Oscar.hasrelshp(H[1], RR) !== nothing
             for g in gens(H[1])
               @test degree(H[2](g)) == grp_elem
               @test (H[2].g)(RR(g)) == g
             end
             @test dim(H[1]) == dim_test #
	  end
        end
        #H_quo = homogeneous_component(R_quo, grp_elem)
        #Oscar.hasrelshp(H_quo[1], R_quo) !== nothing
        #for g in gens(H_quo[1])
        #    degree(H_quo[2](g)) == grp_elem
        #    (H_quo[2].g)(R_quo(g)) == g
        #end
    end
  end
end

@testset "Coercion" begin
  R, (x, y) = grade(PolynomialRing(QQ, ["x", "y"])[1]);
  Q = quo(R, ideal([x^2, y]))[1];
  @test parent(Q(x)) === Q
  @test parent(Q(gens(R.R)[1])) === Q
end

@testset "Evaluation" begin
  R, (x,y) = grade(PolynomialRing(QQ, ["x", "y"])[1]);
  @test x(y, x) == y
end

@testset "Promotion" begin
  R, (x,y) = grade(PolynomialRing(QQ, ["x", "y"])[1]);
  @test x + QQ(1//2) == x + 1//2
end

@testset "Degree" begin
  R, (x,y) = grade(PolynomialRing(QQ, ["x", "y"])[1]);
  @test_throws ArgumentError degree(zero(R))

  Z = abelian_group(0)
  R, (x,y) = filtrate(PolynomialRing(QQ, ["x", "y"])[1], [Z[1], Z[1]], (x,y) -> x[1] > y[1])
  @test degree(x+y^3) == 1*Z[1]
end

@testset "Grading" begin
  R, (x,y) = grade(PolynomialRing(QQ, ["x", "y"])[1]);
  D = grading_group(R)
  @test isisomorphic(D, abelian_group([0]))
end

# Conversion bug

begin
  R, (x,y) = QQ["x", "y"]
  R_ext, _ = PolynomialRing(R, ["u", "v"])
  S, (u,v) = grade(R_ext, [1,1])
  @test S(R_ext(x)) == @inferred S(x)
end

begin
  R, (x,y) = QQ["x", "y"]
  P, (s, t) = R["s", "t"]
  S, (u, v) = grade(P, [1,1])
  @test one(QQ)*u == u
end

begin
  R, (x, y) = QQ["x", "y"]
  S, (u, v) = grade(R)
  @test hash(S(u)) == hash(S(u))

  D = Dict(u => 1)
  @test haskey(D, u)
  @test !haskey(D, v)
end
