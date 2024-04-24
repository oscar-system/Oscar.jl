#cyclic7

using Oscar

using GroebnerWalk
R, (z0, z1, z2, z3, z4, z5, z6) = polynomial_ring(QQ, ["z$index" for index in  0:6 ])

I = ideal([z0 + z1 + z2 + z3 + z4 + z5 + z6,

 z0*z1 + z1*z2 + z2*z3 + z3*z4 + z4*z5 + z5*z6 + z6*z0,

 z0*z1*z2 + z1*z2*z3 + z2*z3*z4 + z3*z4*z5 + z4*z5*z6 + z5*z6*z0 + z6*z0*z1,

 z0*z1*z2*z3 + z1*z2*z3*z4 + z2*z3*z4*z5 + z3*z4*z5*z6 + z4*z5*z6*z0
+ z5*z6*z0*z1 + z6*z0*z1*z2,

 z0*z1*z2*z3*z4 + z1*z2*z3*z4*z5 + z2*z3*z4*z5*z6 + z3*z4*z5*z6*z0 
+ z4*z5*z6*z0*z1 + z5*z6*z0*z1*z2 + z6*z0*z1*z2*z3,

 z0*z1*z2*z3*z4*z5 + z1*z2*z3*z4*z5*z6 + z2*z3*z4*z5*z6*z0 + z3*z4*z5*z6*z0*z1
+ z4*z5*z6*z0*z1*z2 + z5*z6*z0*z1*z2*z3 + z6*z0*z1*z2*z3*z4,

 z0*z1*z2*z3*z4*z5*z6 - 1])

 Ginit = groebner_basis(I)

 set_verbosity_level(:groebner_walk, 1)

t_b = @elapsed Gb = groebner_basis(I, ordering = lex(R), complete_reduction = true) 

t_s = @elapsed Gs = groebner_walk(I, lex(R), algorithm =:standard) 

t_g = @elapsed Gg = groebner_walk(I, lex(R), algorithm =:generic) 

t_p = @elapsed Gp = groebner_walk(I, lex(R), algorithm =:perturbed) 

t_f = @elapsed Gf = groebner_basis(I, ordering = lex(R), algorithm =:fglm) 