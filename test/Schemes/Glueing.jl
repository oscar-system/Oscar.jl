@testset "glueings" begin
  R, (x,y,z) = QQ["x", "y", "z"]
  A3 = Spec(R)
  set_name!(A3, "𝔸³")
  f = (x*y-z^2)
  #f = (x*y-z^2)*(x+y+2*z)
  X = subscheme(A3, f)
  set_name!(X, "X")
  U = SpecOpen(A3, [x,y,z])
  UX = intersect(X, U)
  d = find_non_zero_divisor(UX)
  S, (u,v) = QQ["u", "v"]
  A2 = Spec(S)
  set_name!(A2, "𝔸²")
  f = maximal_extension(X, A2, [x, z//y])
  a = generic_fractions(f)
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

  gxy = Glueing(Ax, Ay, restrict(fxy, domain(fxy), domain(fyx)), restrict(fyx, domain(fyx), domain(fxy)))
  gxz = Glueing(Ax, Az, restrict(fxz, domain(fxz), domain(fzx)), restrict(fzx, domain(fzx), domain(fxz)))
  gyz = Glueing(Ay, Az, restrict(fyz, domain(fyz), domain(fzy)), restrict(fzy, domain(fzy), domain(fyz)))

  gyz_alt = compose(gxy, gxz)
  @test gyz == maximal_extension(gyz_alt)
end

@testset "further glueings" begin
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

  f = SpecMor(Ux, Vu, pbf)
  g = SpecMor(Vu, Ux, pbg)

  simpleG = SimpleGlueing(X, Y, f, g)
  G1 = Glueing(simpleG)
  @test sprint(show, G1) isa String

  Vv = PrincipalOpenSubset(Y, v)
  Wb = PrincipalOpenSubset(Z, b)

  f = SpecMor(Vv, Wb, [u, 1//v])
  g = SpecMor(Wb, Vv, [a, 1//b])
  Vvo = SpecOpen(Vv)
  Wbo = SpecOpen(Wb)
  simpleG2 = SimpleGlueing(Y, Z, f, g)
  @test compose(simpleG, simpleG2) == compose(simpleG, inverse(simpleG2))
  @test compose(inverse(simpleG), simpleG2) == compose(simpleG, inverse(simpleG2))
  @test compose(inverse(simpleG), inverse(simpleG2)) == compose(simpleG, inverse(simpleG2))
  G2 = Glueing(Y, Z, 
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
  @attributes mutable struct DummyGlueing{
                                          LeftSpecType<:AbsSpec, 
                                          RightSpecType<:AbsSpec,
                                          LeftOpenType<:SpecOpen, 
                                          RightOpenType<:SpecOpen,
                                          LeftMorType<:SpecOpenMor,
                                          RightMorType<:SpecOpenMor
                                         } <: AbsGlueing{
                                                         LeftSpecType,
                                                         RightSpecType,
                                                         LeftOpenType,
                                                         RightOpenType,
                                                         LeftMorType,
                                                         RightMorType
                                                        }
    G::Glueing
    function DummyGlueing(G::Glueing{A, B, C, D, E, F}) where {A, B, C, D, E, F}
      return new{A, B, C, D, E, F}(G)
    end
  end

  function Oscar.underlying_glueing(DG::DummyGlueing)
    return DG.G
  end

  function (DG::DummyGlueing)(G::Glueing)
    return DummyGlueing(G)
  end

  ### now everything should work
  DG = DummyGlueing(G1)
  @test (X, Y) == patches(DG)
  @test glueing_morphisms(G1) == glueing_morphisms(DG)
  @test glueing_domains(G1) == glueing_domains(DG)
  @test inverse(DG) == DG
end
