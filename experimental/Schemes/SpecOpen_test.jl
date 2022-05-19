R, (x,y,z) = QQ["x", "y", "z"]

X = Spec(R)

N = subscheme(X, [x,y])

U = complement(X, N)

V = SpecOpen(X, [x^2, y-x])

I = ideal(R, [x*y-z^2])

Y = subscheme(X, I)

U = SpecOpen(Y, [x^5, y-7*x^7])

f = maximal_extension(Y, x//z)

S, (u, v) = PolynomialRing(QQ, ["u", "v"])

Z = Spec(S)

g = maximal_extension(Y, Z, [x//z, y])

g = restriction(g, domain(g), SpecOpen(Z, [v]))
