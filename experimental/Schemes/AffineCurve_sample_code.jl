R, (x,y) = QQ["x", "y"]
C = AffineCurve(R, [x], [y])
A = Spec(R)
pb = hom(OO(A), OO(C), gens(OO(C)))
f = SpecMor(C, A, pb)
