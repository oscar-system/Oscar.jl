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
