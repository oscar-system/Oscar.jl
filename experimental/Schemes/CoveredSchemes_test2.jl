R, (x,y,z,w) = QQ["x", "y", "z", "w"]
Rx, (yx,zx,wx) = QQ["y/x", "z/x", "w/x"]
Ry, (xy,zy,wy) = QQ["x/y", "z/y", "w/y"]
Rz, (xz,yz,wz) = QQ["x/z", "y/z", "w/z"]
Rw, (xw,yw,zw) = QQ["x/w", "y/w", "z/w"]

px = AlgebraHomomorphism(R, Rx, [one(Rx), yx, zx, wx])
py = AlgebraHomomorphism(R, Ry, [xy, one(Ry), zy, wy])
pz = AlgebraHomomorphism(R, Rz, [xz, yz, one(Rz), wz])
pw = AlgebraHomomorphism(R, Rw, [xw, yw, zw, one(Rw)])

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

gxy = SpecMor(hypersurface_complement(UXx, yx), hypersurface_complement(UXy, xy), [1//yx, zx//yx, wx//yx])
Gxy = Glueing(UXx, UXy, gxy)
gxz = SpecMor(hypersurface_complement(UXx, zx), hypersurface_complement(UXz, xz), [1//zx, yx//zx, wx//zx])
Gxz = Glueing(UXx, UXz, gxz)
gxw = SpecMor(hypersurface_complement(UXx, wx), hypersurface_complement(UXw, xw), [1//wx, yx//wx, zx//wx])
Gxw = Glueing(UXx, UXw, gxw)

Gyxz = compose(Gxy, Gxz)
Gyz_by_ext = maximal_extension(Gyxz)[1]

gyz = SpecMor(hypersurface_complement(UXy, zy), hypersurface_complement(UXz, yz), [xy//zy, 1//zy, wy//zy])
Gyz = Glueing(UXy, UXz, gyz)

Gyz == Gyz_by_ext
