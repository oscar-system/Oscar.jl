
using Oscar

R, (x,y,z,t,u,v) = QQ[:x,:y,:z,:t,:u,:v]

I = ideal([
    2*x^2+2*y^2+2*z^2+2*t^2+2*u^2+v^2-v,
    x*y+y*z+2*z*t+2*t*u+2*u*v-u,
    2*x*z+2*y*t+2*z*u+u^2+2*t*v-t,
    2*x*t+2*y*u+2*t*u+2*z*v-z,
    t^2+2*x*v+2*y*v+2*z*v-y,
    2*x+2*y+2*z+2*t+2*u+v-1
])

o_s = degrevlex(R)
o_t = lex(R)

set_verbosity_level(:groebner_walk, 0)

t_s = @elapsed Gs = groebner_walk(I, lex(R), algorithm=:standard)
t_g = @elapsed Gg = groebner_walk(I, lex(R), algorithm=:generic)
t_b = @elapsed Gb = groebner_basis(I, ordering=lex(R))

