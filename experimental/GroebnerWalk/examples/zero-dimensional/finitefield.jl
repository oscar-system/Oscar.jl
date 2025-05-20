#=
  Example over finite field, from `Groebner.jl` docstring line 1232
=#
using Oscar
R, (x1, x2, x3, x4) = polynomial_ring(GF(32003), [:x1, :x2, :x3, :x4])

J = ideal(R, [
  x1 + 2 * x2 + 2 * x3 + 2 * x4 - 1,
  x1^2 + 2 * x2^2 + 2 * x3^2 + 2 * x4^2 - x1,
  2 * x1 * x2 + 2 * x2 * x3 + 2 * x3 * x4 - x2,
  x2^2 + 2 * x1 * x3 + 2 * x2 * x4 - x3
])

t_s = @elapsed Gs = groebner_walk(J, lex(R), algorithm=:standard) #4.11
t_g = @elapsed Gg = groebner_walk(J, lex(R), algorithm=:generic) #0.8s
