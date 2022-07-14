# S is NS of a K3 surface with finite automorphism discriminant_group
# It has exactly 4 (-2)-curves
# and the roots are the standard basis vectors.
# Draw the dual graph
# It has one symmetry. It is of order 2.
add_verbose_scope(:K3Auto)
set_verbose_level(:K3Auto,5)

gram = QQ[-2 1 0 0; 1 -2 1 1; 0 1 -2 1; 0 1 1 -2]
S = Zlattice(gram=gram)
# Construct an embedding of S into L \cong L_{10}.
R = rescale(root_lattice(:E,6),-1)
SR,iS,iR = orthogonal_sum(S,R)
V = ambient_space(SR)
S = lattice(V,basis_matrix(S)*iS.matrix)
R = lattice(V,basis_matrix(R)*iR.matrix)
g = gens(discriminant_group(SR))
i = g[1]+g[2]
#@assert Hecke.quadratic_product(i)==0 #bug
Bi = vcat(basis_matrix(SR),matrix(QQ,1,10,lift(i)))
L = lattice(V,Bi,isbasis=false)

# Find an isomorphism L \cong U + E_8
f = QQ[0 1 1 1 0 0 0 0 0 0]
z = QQ[1 0 0 0 0 0 0 0 0 0]
@assert inner_product(V,f,f)==0
U = lattice(V,vcat(f,z))
E8 = Hecke.orthogonal_submodule(L,U)
e8 = rescale(root_lattice(:E,8),-1)

# normalize the basis
_,T = isisometric(e8,E8, ambient_representation=false)
E8 = lattice(ambient_space(E8),T*basis_matrix(E8))
@assert gram_matrix(E8) == gram_matrix(e8)
B = vcat(f, z, basis_matrix(E8))
UE8 = lattice(ambient_space(E8),B)
Bdual = inv(gram_matrix(ambient_space(UE8))*transpose(B))

# this one does not have ample projection
weyl = QQ[30 1 1 1 1 1 1 1 1 1]*Bdual
h = ZZ[3 12 11 11]*basis_matrix(S)  #an ample vector
# confirm that h is in the interior of a weyl chamber,
# i.e. check that Q does not contain any -2 vector and h^2>0
Q = Hecke.orthogonal_submodule(S, lattice(V, h))
@test minimum(rescale(Q, -1)) > 2
# make weylS ample
weylS = hcat(weyl[1,1:4],zero_matrix(QQ,1,6))
separating_roots = oscar.alg23(S, h, weylS, -2)
# order the reflections appropriately. I did not check the math behind this. So it might be wrong.
sort!(separating_roots, by=r->inner_product(V,r,h)[1,1]//inner_product(V,r,weyl)[1,1])
for r in separating_roots
  global weyl = weyl + inner_product(V, weyl, r)*r
end
weylS = hcat(weyl[1,1:4],zero_matrix(QQ,1,6))
@test length(oscar.alg23(S, h, weylS, -2))==0
@assert weyl == QQ[30 61 146//3 136//3 11//3 -11//3 3 -7//3 -20//3 -1]
# the weyl vector is S-degenerate
# make weyl S-nondegenerate as follows
weyl,_,_ = oscar.nondeg_weyl_new(L,S,weyl,weyl,h)
# since the output is not deterministic we hardcode:
#weyl = QQ[30   61   51   42   5   2   -4   -5   -3   -1]

Gamma, W, DD, B = oscar.alg61(L,S,weyl)


# Another example with finite automorphism group
S,iU,_=orthogonal_sum(Zlattice(gram=ZZ[0 1; 1 -2]),rescale(root_lattice(:D,4),-1))
L,S,iS, R,iR = oscar.embed_in_unimodular(S::ZLat, 10)
V = ambient_space(L)
# find a hyperbolic plane
U = lattice(V, iU.matrix*iS.matrix)
weyl,_ = oscar.weyl_vector(L, U)

h = ZZ[25 6 -10 -9 -16 -10]*basis_matrix(S)  #an ample vector
@assert inner_product(V,h,h)[1,1]>0
@assert all([a>0 for a in inner_product(V,h,basis_matrix(S))])
# confirm that h is in the interior of a weyl chamber,
# i.e. check that Q does not contain any -2 vector and h^2>0
Q = Hecke.orthogonal_submodule(S, lattice(V, h))
@test minimum(rescale(Q, -1)) > 2

weyl,_,_ = oscar.nondeg_weyl_new(L,S,weyl,weyl,h)

Gamma, W, DD, B = oscar.alg61(L,S,weyl)

# another example
S,iU,_=orthogonal_sum(Zlattice(gram=ZZ[0 1; 1 -2]),Zlattice(gram=ZZ[-50;]))
L,S,iS, R,iR = oscar.embed_in_unimodular(S::ZLat, 18)
V = ambient_space(L)
# find a hyperbolic plane
U = lattice(V, iU.matrix*iS.matrix)
weyl,u0 = oscar.weyl_vector(L, U)

h = ZZ[40 1 -1 ]*basis_matrix(S)  #an ample vector
@assert inner_product(V,h,h)[1,1]>0
@assert all([a>0 for a in inner_product(V,h,basis_matrix(S))])
# confirm that h is in the interior of a weyl chamber,
# i.e. check that Q does not contain any -2 vector and h^2>0
Q = Hecke.orthogonal_submodule(S, lattice(V, h))
@test minimum(rescale(Q, -1)) > 2

weyl,_ = oscar.nondeg_weyl_new(L,S,weyl,weyl,h)

Gamma, W, DD, B = oscar.alg61(L,S,weyl)



C = lattice(V,common_invariant(Gamma)[2])
diagonal(rational_span(C))


zero_entropy_candidates = oscar.parse_zero_entropy()

zero_entropy_candidates = [Zlattice(gram=g) for g in zero_entropy_candidates]

S = zero_entropy_candidates[3]

