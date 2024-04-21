using Oscar 
using GroebnerWalk
R, (x,y,z) = polynomial_ring(QQ, ["x","y","z"])

I = ideal([6 + 3*x^3+16*x^2*z+14*x^2*y^3, 6+y^3*z+17*x^2*z^2+7*x*y^2*z^2+13*x^3*z^2])

# This is an example 3.8 from Tran 2004. 
# In M2, "gb I" takes ~90s
# standard Gröbner Walk takes ~12s
# generic Gröbner Walk takes ~2s 

t_buchberger = @elapsed groebner_basis(I, ordering = lex(R)) 
#this takes longer than 30 mins in OSCAR, compared with ~90s in M2
#Is this to be expected? 


t_standardwalk = @elapsed groebner_walk(I, lex(R); algorithm =:standard) 

t_generic = @elapsed groebner_walk(I, lex(R); algorithm =:generic) #takes more than 1hr...


#the line above terminates in a few secs 

#you should count the number of conversions

#problem: the two functions above appear to do exactly the same thing! 
#"Results for standard_walk....etc. 

