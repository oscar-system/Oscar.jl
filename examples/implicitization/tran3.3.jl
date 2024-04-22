#Example 3.3 from Tran
# This is an example 3.3 from "A fast algorithm for Gröbner basis conversion and its applications"  ( Tran 2000)
# 1-dimensional ideal (application not given)
using Oscar 
using GroebnerWalk
R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

I = ideal([16 + 3*x^3+16*x^2*z+14*x^2*y^3, 6+y^3*z+17*x^2*z^2+7*x*y^2*z^2+13*x^3*z^2])


Gdegrevlex = groebner_basis(I, complete_reduction = true)

set_verbosity_level(:groebner_walk, 1)

t_b = @elapsed Gb = groebner_basis(I, ordering = lex(R), complete_reduction = true) #doesn't terminate after 90mins

t_s = @elapsed Gs = groebner_walk(I, lex(R), algorithm =:standard) #~10s faster than m2! (but remember; no reduction! ))

t_g = @elapsed Gg = groebner_walk(I, lex(R), algorithm =:generic) #87s. why so long? we have many conversions, maybe inflating weight vectors? 

t_p = @elapsed Gp = groebner_walk(I, lex(R), algorithm =:perturbed) #54s, long final conversion. Also, I don't understand the output
