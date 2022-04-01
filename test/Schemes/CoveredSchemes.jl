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

  gxy = Glueing(Ax, Ay, restriction(fxy, domain(fxy), domain(fyx)), restriction(fyx, domain(fyx), domain(fxy)))
  gxz = Glueing(Ax, Az, restriction(fxz, domain(fxz), domain(fzx)), restriction(fzx, domain(fzx), domain(fxz)))
  gyz = Glueing(Ay, Az, restriction(fyz, domain(fyz), domain(fzy)), restriction(fzy, domain(fzy), domain(fyz)))

  gyz_alt = compose(gxy, gxz)
  @test gyz == maximal_extension(gyz_alt)
end

@testset "covered schemes 1" begin
  # manual glueing of ℙ³ and a quadric surface in there
  R, (x,y,z,w) = QQ["x", "y", "z", "w"]
  Rx, (yx,zx,wx) = QQ["y/x", "z/x", "w/x"]
  Ry, (xy,zy,wy) = QQ["x/y", "z/y", "w/y"]
  Rz, (xz,yz,wz) = QQ["x/z", "y/z", "w/z"]
  Rw, (xw,yw,zw) = QQ["x/w", "y/w", "z/w"]

  px = hom(R, Rx, [one(Rx), yx, zx, wx])
  py = hom(R, Ry, [xy, one(Ry), zy, wy])
  pz = hom(R, Rz, [xz, yz, one(Rz), wz])
  pw = hom(R, Rw, [xw, yw, zw, one(Rw)])

  F = x*w-y*z
  I = ideal(R, F)
  IP = Spec(R)
  IPX = subscheme(IP, I)
  Ux = Spec(Rx)
  UXx = subscheme(Ux, ideal(Rx, px(F)))
  Uy = Spec(Ry)
  UXy = subscheme(Uy, ideal(Ry, py(F)))
  Uz = Spec(Rz)
  UXz = subscheme(Uz, ideal(Rz, pz(F)))
  Uw = Spec(Rw)
  UXw = subscheme(Uw, ideal(Rw, pw(F)))

  Gxy = Glueing(UXx, yx, UXy, xy, [1//yx, zx//yx, wx//yx], [1//xy, zy//xy, wy//xy])

  Gxz = Glueing(UXx, zx, UXz, xz, [1//zx, yx//zx, wx//zx], [yz//xz, 1//xz, wz//xz])
  Gxw = Glueing(UXx, wx, UXw, xw, [1//wx, yx//wx, zx//wx], [yw//xw, zw//xw, 1//xw])

  Gyxz = compose(Gxy, Gxz)
  Gyz_by_ext = maximal_extension(Gyxz)
  @test Gyxz == Gyz_by_ext
end
