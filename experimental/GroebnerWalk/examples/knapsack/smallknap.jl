#smallknap

using Oscar 
R, (t, x1, x2, x3) = polynomial_ring(QQ, [:t,:x1, :x2, :x3])

f1 = x1 - t^5
f2 = x2 - t^12
f3 = x3 - t^20

I = ideal([f1, f2, f3])

set_verbosity_level(:groebner_walk, 1)

o_t = weight_ordering([1,0,0,0], degrevlex(R))
o_s = weight_ordering([0,1,1,1], degrevlex(R))

Gg = groebner_walk(I, o_t, o_s, algorithm =:generic)

