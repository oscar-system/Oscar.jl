@testset "gluings" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  A3 = Spec(R)
  set_name!(A3, "𝔸³")
  f = (x*y-z^2)
  #f = (x*y-z^2)*(x+y+2*z)
  X = subscheme(A3, f)
  set_name!(X, "X")
  U = SpecOpen(A3, [x,y,z])
  UX = intersect(X, U)
  d = Oscar.find_non_zero_divisor(UX)
  S, (u,v) = QQ["u", "v"]
  A2 = Spec(S)
  set_name!(A2, "𝔸²")
  f = maximal_extension(X, A2, [x, z//y])
  a = Oscar.generic_fractions(f)
  @test maximal_extension(X, A2, a) == f

  Sx, (yx, zx) = QQ["yx", "zx"]
  Sy, (xy, zy) = QQ["xy", "zy"]
  Sz, (xz, yz) = QQ["xz", "yz"]

  Ax = Spec(Sx)
  Ay = Spec(Sy)
  Az = Spec(Sz)

  fxy = maximal_extension(Ax, Ay, [1//yx, zx//yx])
  fyx = maximal_extension(Ay, Ax, [1//xy, zy//xy])
  fxz = maximal_extension(Ax, Az, [1//zx, yx//zx])
  fzx = maximal_extension(Az, Ax, [yz//xz, 1//xz])
  fyz = maximal_extension(Ay, Az, [xy//zy, 1//zy])
  fzy = maximal_extension(Az, Ay, [xz//yz, 1//yz])

  gxy = Gluing(Ax, Ay, restrict(fxy, domain(fxy), domain(fyx)), restrict(fyx, domain(fyx), domain(fxy)))
  gxz = Gluing(Ax, Az, restrict(fxz, domain(fxz), domain(fzx)), restrict(fzx, domain(fzx), domain(fxz)))
  gyz = Gluing(Ay, Az, restrict(fyz, domain(fyz), domain(fzy)), restrict(fzy, domain(fzy), domain(fyz)))

  gyz_alt = compose(gxy, gxz)
  @test gyz == maximal_extension(gyz_alt)
end

@testset "further gluings" begin
  R, (x, y) = QQ["x", "y"]
  S, (u, v) = QQ["u", "v"]
  T, (a, b) = QQ["a", "b"]

  X = Spec(R)
  Y = Spec(S)
  Z = Spec(T)

  Ux = PrincipalOpenSubset(X, x)
  Vu = PrincipalOpenSubset(Y, u)

  pbf = hom(OO(Vu), OO(Ux), [inv(OO(Ux)(x)), OO(Ux)(y)])
  pbg = hom(OO(Ux), OO(Vu), [inv(OO(Vu)(u)), OO(Vu)(v)])

  f = morphism(Ux, Vu, pbf)
  g = morphism(Vu, Ux, pbg)

  simpleG = SimpleGluing(X, Y, f, g)
  G1 = Gluing(simpleG)
  @test sprint(show, G1) isa String

  Vv = PrincipalOpenSubset(Y, v)
  Wb = PrincipalOpenSubset(Z, b)

  f = morphism(Vv, Wb, [u, 1//v])
  g = morphism(Wb, Vv, [a, 1//b])
  Vvo = SpecOpen(Vv)
  Wbo = SpecOpen(Wb)
  simpleG2 = SimpleGluing(Y, Z, f, g)
  @test compose(simpleG, simpleG2) == compose(simpleG, inverse(simpleG2))
  @test compose(inverse(simpleG), simpleG2) == compose(simpleG, inverse(simpleG2))
  @test compose(inverse(simpleG), inverse(simpleG2)) == compose(simpleG, inverse(simpleG2))
  G2 = Gluing(Y, Z, 
               SpecOpenMor(Vvo, Wbo, [compose(f, inclusion_morphism(Wb, Z))]), 
               SpecOpenMor(Wbo, Vvo, [compose(g, inclusion_morphism(Vv, Y))]))

  G3 = compose(G1, G2)
  @test G3 == compose(G1, inverse(G2))
  @test G3 == inverse(compose(inverse(G1), G2))
  @test G3 == inverse(G3)

  Xsub = subscheme(X, y-x^2)
  Ysub = subscheme(Y, u^2*v-1)
  G1res = restrict(G1, Xsub, Ysub)

  ### test the abstract interface
  @attributes mutable struct DummyGluing{
                                          LeftSpecType<:AbsSpec, 
                                          RightSpecType<:AbsSpec,
                                          LeftOpenType<:SpecOpen, 
                                          RightOpenType<:SpecOpen,
                                          LeftMorType<:SpecOpenMor,
                                          RightMorType<:SpecOpenMor
                                         } <: AbsGluing{
                                                         LeftSpecType,
                                                         RightSpecType,
                                                         LeftOpenType,
                                                         RightOpenType,
                                                         LeftMorType,
                                                         RightMorType
                                                        }
    G::Gluing
    function DummyGluing(G::Gluing{A, B, C, D, E, F}) where {A, B, C, D, E, F}
      return new{A, B, C, D, E, F}(G)
    end
  end

  function Oscar.underlying_gluing(DG::DummyGluing)
    return DG.G
  end

  function (DG::DummyGluing)(G::Gluing)
    return DummyGluing(G)
  end

  ### now everything should work
  DG = DummyGluing(G1)
  @test (X, Y) == patches(DG)
  @test gluing_morphisms(G1) == gluing_morphisms(DG)
  @test gluing_domains(G1) == gluing_domains(DG)
  @test inverse(DG) == DG
end

@testset "base change" begin
  kk, pr = quo(ZZ, 5)
  IP1 = covered_scheme(projective_space(ZZ, 1))
  C = default_covering(IP1)
  G = first(values(gluings(C)))
  GG = base_change(pr, G)
  @test underlying_gluing(GG) isa SimpleGluing
end

@testset "extra lazy gluing domains" begin
  P3 = projective_space(QQ, 3)
  S = homogeneous_coordinate_ring(P3)
  (x, y, z, w) = gens(S)
  I = ideal(S, [x*z - y*w, x^4 + y^4 + z^4 + w^4])
  Y = covered_scheme(P3)
  II = ideal_sheaf(P3, I)
  pr = blow_up(II)
  X = domain(pr)

  for G in values(gluings(Oscar.simplified_covering(X)))
    @test !isdefined(G, :G)
    gluing_domains(G)
    @test !isdefined(G, :G)
    @test isdefined(G, :gluing_domains)
  end

  for G in values(gluings(Oscar.simplified_covering(X)))
    @test !isdefined(G, :G)
    underlying_gluing(G)
    @test isdefined(G, :G)
    @test isdefined(G, :gluing_domains)
    @test (G.gluing_domains[1] === gluing_domains(underlying_gluing(G))[1])
    @test (G.gluing_domains[2] === gluing_domains(underlying_gluing(G))[2])
  end
end

