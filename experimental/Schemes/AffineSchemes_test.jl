R, rvars = QQ["x", "y"]
x = rvars[1]
y = rvars[2]

S, svars = QQ["v", "w"]
v = svars[1]
w = svars[2]

X = Spec(R)
U = hypersurface_complement(X, x)
Y = Spec(S)
V = hypersurface_complement(Y, v)
P = subscheme(Y, v*w)

#imgs = localized_ring(OO(U)).([1//x, y])
imgs = OO(U).([1//x, y])
@show imgs
phi = MPolyQuoLocalizedRingHom(OO(V), OO(U), [1//x, y])

f = SpecHom(U, V, phi)
preimage(f, P)
