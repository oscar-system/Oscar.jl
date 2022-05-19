R, (x,y) = QQ["x", "y"]
S, (u,v) = QQ["u", "v"] 
T, (a,b) = QQ["a", "b"]
I = ideal(R, 2*x^2+3*y^2-4)
X = Spec(R)
set_name!(X, "X")
Y = Spec(S)
set_name!(Y, "Y")
Z = Spec(T)
set_name!(Z, "Z")
U = hypersurface_complement(X, y)
set_name!(U, "U")
V1 = hypersurface_complement(Y, v)
set_name!(V1, "V1")
V2 = hypersurface_complement(Y, u)
set_name!(V2, "V2")
W = hypersurface_complement(Z, a)
set_name!(W, "W")

ZX = subscheme(X, I)
ZU = intersect(ZX, U)

f = SpecMor(U, V1, [x//y, 1//y])
gamma = Glueing(X, Y, f)

g = SpecMor(V2, W, [1//u, v//u])
delta = Glueing(Y, Z, g)

epsilon = compose(gamma, delta)

h = maximal_extension(epsilon)[1]

CX = subscheme(X, I)
CU = intersect(CX, U)
CV1 = preimage(inverse(f), CU)
CY = closure(CV1, Y)
CV2 = intersect(CY, V2)
CW = preimage(inverse(g), CV2)
CZ = closure(CW, Z)

CU2 = intersect(CX, domain(glueing_morphism(h)))
CW2 = preimage(inverse(glueing_morphism(h)), CU2)
CZ2 = closure(CW2, Z)

CZ == CZ2


