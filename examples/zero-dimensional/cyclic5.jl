using Oscar
using GroebnerWalk

R, (a, b, c, d, x)  = polynomial_ring(QQ, ["a", "b", "c", "d", "x"])

I = ideal([ a + b + c + d + x,
a*b + b*c + c*d + d*x + x*a,
a*b*c + b*c*d + c*d*x + d*x*a + x*a*b,
a*b*c*d + b*c*d*x + c*d*x*a + d*x*a*b + x*a*b*c,
a*b*c*d*x - 1] )

#G = groebner_basis(I, ordering = lex(R), complete_reduction = true)
o_s = degrevlex(R)
o_t = lex(R)
G1 = groebner_basis(I, ordering = degrevlex(R), complete_reduction = true)


set_verbosity_level(:groebner_walk, 1)

t_b = @elapsed Gb = groebner_basis(I, ordering = lex(R), complete_reduction = true) 

t_s = @elapsed Gs = groebner_walk(I, lex(R), algorithm =:standard) 

t_g = @elapsed Gg = groebner_walk(I, lex(R), algorithm =:generic) 

t_p = @elapsed Gp = groebner_walk(I, lex(R), algorithm =:perturbed) 

t_f = @elapsed Gf = groebner_basis(I, ordering = lex(R), algorithm =:fglm) 
#= What do we learn from this example? 

Gg1. standardwalk is fastest, but it doesn't give us a REDUCED Gröbner basis
2. the output of genericwalk and buchberger (with complete_reduction) is the same UP TO SCALING
3. The elements of the genericwalkbasis have rational coefficients (due to divides()...), whereas Buchberger's output doesnt

=# 