using Oscar

R, (a, b, c, d, x) = polynomial_ring(QQ, [:a, :b, :c, :d, :x])

I = ideal([
  a + b + c + d + x,
  a * b + b * c + c * d + d * x + x * a,
  a * b * c + b * c * d + c * d * x + d * x * a + x * a * b,
  a * b * c * d + b * c * d * x + c * d * x * a + d * x * a * b + x * a * b * c,
  a * b * c * d * x - 1
])

o_s = degrevlex(R)
o_t = lex(R)

set_verbosity_level(:groebner_walk, 1)

Gs = groebner_walk(I, lex(R), algorithm=:standard)
Gg = groebner_walk(I, lex(R), algorithm=:generic)

