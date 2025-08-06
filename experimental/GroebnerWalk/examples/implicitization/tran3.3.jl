#= 
   Example 3.3 from Tran. "A fast algorithm for Gröbner basis conversion and its applications" (2000)
   A 1-dimensional ideal (application not given)
=#
using Oscar 
R, (x,y,z) = polynomial_ring(QQ, [:x,:y,:z])
I = ideal([16 + 3*x^3+16*x^2*z+14*x^2*y^3, 6+y^3*z+17*x^2*z^2+7*x*y^2*z^2+13*x^3*z^2])

set_verbosity_level(:groebner_walk, 1)
Gs = groebner_walk(I, lex(R), algorithm =:standard)

