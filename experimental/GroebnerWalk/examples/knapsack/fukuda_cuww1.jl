#=
  Hard integer knapsack problem from Aardal, Karen and Lenstra. ‘Hard Equality Constrained Integer Knapsacks’. (2004)
=#

using Oscar
R, (t, x1, x2, x3, x4, x5) = polynomial_ring(QQ, [:t, :x1, :x2, :x3, :x4, :x5])

f = 12223 * x1 + 12224 * x2 + 36674 * x3 + 61119 * x4 + 85569 * x5 - 89643481

f1 = x1 - t^1223
f2 = x2 - t^1224
f3 = x3 - t^36674
f4 = x4 - t^61119
f5 = x5 - t^85569
I = ideal([f1, f2, f3, f4, f5])

set_verbosity_level(:groebner_walk, 1)

o_t = weight_ordering([1, 0, 0, 0, 0, 0], degrevlex(R))
o_s = weight_ordering([0, 1, 1, 1, 1, 1], degrevlex(R))

ts = @elapsed Gs = groebner_walk(I, o_t, o_s, algorithm=:standard)

