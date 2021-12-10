R, (x,y) = QQ["x", "y"]

M0 = free_module(R, 1)
M1 = free_module(R, 1)

X = hom(M0, M1, [x*gen(M1, 1)])
Y = hom(M0, M1, [y*gen(M1, 1)])

C = BoundedCocomplex([M0, M1], [X], shift=2)
f = BoundedCocomplexHom(C, C, [hom(M, M, [y*e for e in gens(M)]) for M in cochains(C)])
(Cf, i, p) = mapping_cone(f)
