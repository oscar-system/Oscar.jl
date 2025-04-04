# definition of monoid algebra as quotient of polynomial ring
S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
J = ideal(S,[x*z-y^2])
R_Q,phi = quo(S,J)
x,y,z = gens(R_Q)

# get MonoidAlgebra
kQ = get_monoid_algebra(R_Q)

#example of computing local cohomology modules 
I_M = ideal(kQ,[x^2*z,x^4*y])
I = ideal(kQ,[x,y,z]) #maximal ideal
M = quotient_ring_as_module(I_M)

# compute local cohomology modules supported on I = (x,y,z)
H0 = zeroth_local_cohomology(quotient_ring_as_module(I_M),I)
H1 = local_cohomology(I_M,I,1)
H2 = local_cohomology(I_M,I,2)
H3 = local_cohomology(I_M,I,3)
H4 = local_cohomology(I_M,I,4)

# two sectors of local_cohomology
H1.sectors[1]
H1.sectors[20]

# check if local cohomology is zero
is_zero(H0)
lc_zero(H1)
lc_zero(H2)
lc_zero(H3)
lc_zero(H4)