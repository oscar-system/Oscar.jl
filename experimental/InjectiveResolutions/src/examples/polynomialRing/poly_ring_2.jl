# definition of polynomial ring k[x,y]
R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])

FM = graded_free_module(R_Q,2)
ARL = R_Q[x 0; 0 y]
BRL = R_Q[y^4 0; x^2*y -y^3; x^3 -x*y^2; 0 x^4]

_M = SubquoModule(FM,ARL,BRL)

F = graded_free_module(R_Q,1)
# F = free_module(R_Q,1)
Fm = quo_object(F,[x*F[1],y*F[1]])

hom(Fm,_M)
hom(Fm,_M,algorithm=:matrices)

# get MonoidAlgebra
kQ = get_monoid_algebra(R_Q)

M = get_monoid_algebra_module(kQ,_M)

irreducible_res(M) # error...˛


### hom issues
R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])

# using graded_free_module
FM = graded_free_module(R_Q,2)
ARL = R_Q[x 0; 0 y]
BRL = R_Q[y^4 0; x^2*y -y^3; x^3 -x*y^2; 0 x^4]
_M = SubquoModule(FM,ARL,BRL)

F = graded_free_module(R_Q,1)
Fm = quo_object(F,[x*F[1],y*F[1]])

H = hom(Fm,_M) # not working...
H = hom(Fm,_M,algorithm=:matrices)

# using free_module
FM = free_module(R_Q,2)
ARL = R_Q[x 0; 0 y]
BRL = R_Q[y^4 0; x^2*y -y^3; x^3 -x*y^2; 0 x^4]
_M = SubquoModule(FM,ARL,BRL)

F = free_module(R_Q,1)
Fm = quo_object(F,[x*F[1],y*F[1]])
H = hom(Fm,_M)
H = hom(Fm,_M,algorithm=:matrices)