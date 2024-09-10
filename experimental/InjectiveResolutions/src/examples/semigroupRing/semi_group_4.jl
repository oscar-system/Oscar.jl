#---------------------------------------
#------------ preliminary definitions


# definition of semigroup ring
S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
J = ideal(S,[x*z-y^2])
R_Q,phi = quo(S,J)

##define all faces of Q as polyhedron using systems of linear inequalities 
#F0 = cone{(0,0)}
A0 = [1 0;-1 0;0 1;0 -1]
b0 = [0,0,0,0]
F0 = polyhedron(A0,b0)

#F1 = cone{(0,1)}
A1 = [1 0; -1 0; 0 -1]
b1 = [0,0,0]
F1 = polyhedron(A1,b1)

#F2 = cone{(2,1)}
A2 = [1 -2; -1 2; -1 0; 0 -1]
b2 = [0,0,0,0]
F2 = polyhedron(A2,b2)

##define hyperplanes bounding Q as polyhedron
#H_1 hyperplane bounding Q (intersects F1 non-trivially)
A_H1 = [1 0;-1 0]
b_H1 = [0,0]
H1 = polyhedron(A_H1,b_H1)

#H_2 hyperplane bounding Q (intersects F2 non-trivially)
A_H2 = [1 -2; -1 2]
b_H2 = [0,0]
H2 = polyhedron(A_H2,b_H2)

F = [F1,F2] #faces bounding Q (facets)
P_Q = convex_hull(F1,F2)
h1 = FaceQ(ideal(R_Q,[zero(R_Q)]),H1,A_H1,b_H1)
h2 = FaceQ(ideal(R_Q,[zero(R_Q)]),H2,A_H2,b_H2)
H = [h1,h2] #hyperplanes bounding Q

## define zonotope of Q
R1 = convex_hull([0 0;0 1])
R2 = convex_hull([0 0;2 1])
G = R1 + R2 


#-------------------------------------------------
#--------- Example 4
#-------------------------------------------------
I = ideal(R_Q,[x^2*z])
compute_bass_numbers(I,4)
c = [3,3]
P = get_all_ass_primes(I)
m_c = monomial_basis(R_Q,c)[1]
Ic = ideal(R_Q,[m_c])*I
M0,_ = quotient_by_ideal(Ic)
M0c,_ = sub(M0,[m_c*M0[1]])

#compute irreducible resolution
res = irreducible_res(M0c,P,P_Q,G,F,H)

res.irrSums[1].components

res.irrSums[2].components

res.irrSums[3].components

matrix(res.cochainMaps[1])

matrix(res.cochainMaps[2])

matrix(res.cochainMaps[3])

# check if irreducible resolution
length(res.irrSums)
image(res.cochainMaps[1])[1] == kernel(res.cochainMaps[2])[1]
image(res.cochainMaps[2])[1] == kernel(res.cochainMaps[3])[1]
is_injective(res.inclusions[1])
is_injective(res.inclusions[2])
is_injective(res.inclusions[3])
is_surjective(res.inclusions[3])