#---------------------------------------
#------------ preliminary definitions

# definition of polynomial ring k[x,y]
R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])


## faces of Q as polyhedron
# F1 = cone((1,0))
A1 = [-1 0; 1 0;0 -1]
b1 = [0,0,0]
F1 = polyhedron(A1,b1)

# F2 = cone((0,1))
A2 = [-1 0; 0 -1; 0 1]
b2 = [0,0,0]
F2 = polyhedron(A2,b2)

## hyperplanes bounding Q as polyhedron
A_H1 = [-1 0;1 0]
b_H1 = [0,0]
H1 = polyhedron(A_H1,b_H1)

A_H2 = [0 -1;0 1]
b_H2 = [0,0]
H2 = polyhedron(A_H2,b_H2)

F = [F1,F2] #faces bounding Q (facets)
P_Q = convex_hull(F1,F2) # semigroup Q as polyhedron
h1 = FaceQ(ideal(R_Q,[zero(R_Q)]),H1,A_H1,b_H1)
h2 = FaceQ(ideal(R_Q,[zero(R_Q)]),H2,A_H2,b_H2)
H = [h1,h2]

### primitive integer vectors along rays of Q
R1 = convex_hull([0 0;0 1])
R2 = convex_hull([0 0;1 0])
G = R1 + R2 # zonotope as im Lemma 3.10. (HM2004)
c = [1,1] # ZZ^d degree of sum R1 + R2


#-------------------------------------------------
#--------- Example 1 (Example 11.3. in MS2005)
#-------------------------------------------------

I = ideal(R_Q,[x^4,x^2*y^2,y^4])
M,_ = quotient_by_ideal(I)
P = get_all_ass_primes(I) 

#compute an irreducible resolution 
res = irreducible_res(M,P,P_Q,G,F,H)

res.irrSums[1].components

res.irrSums[2].components

matrix(res.cochainMaps[1])

matrix(res.cochainMaps[2])


#check if irreducible resolution
length(res.irrSums)
image(res.cochainMaps[1])[1] == kernel(res.cochainMaps[2])[1]
is_injective(res.inclusions[1])
is_injective(res.inclusions[2])
is_surjective(res.inclusions[2])