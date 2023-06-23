
E8 = rescale(root_lattice(:E, 8),-1)
A2 = rescale(root_lattice(:A, 2),-1)
U = integer_lattice(gram=ZZ[0 1;1 -2])

NS,_ = direct_sum([U,E8,E8,A2])
V = ambient_space(NS)

e = matrix(ZZ,1,20,ones(Int,20))
e[1,1]=51

h = e*inv(gram_matrix(NS))

L, S, weyl = Oscar.borcherds_method_preprocessing(NS, 26, ample=h1)

Xdata, Xaut, Xchambers, Xrational_curves = borcherds_method(L, S, weyl; compute_OR=true)

