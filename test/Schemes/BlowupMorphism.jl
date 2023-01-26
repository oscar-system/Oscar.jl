R, (x,y,z) = QQ["x", "y", "z"]
f = x^2 + y^3 + z^5
X = CoveredScheme(Spec(R, ideal(R, f)))
U = X[1][1] # the first chart

IZ = IdealSheaf(X, U, OO(U).([x, y, z]))

Y = blow_up(IZ)
