# "custom" integer knapsack problem

using Oscar 
R, (t, x1, x2, x3, x4, x5) = polynomial_ring(QQ, [:t,:x1, :x2, :x3, :x4, :x5])






f1 = x1 - t^329221
f2 = x2 - t^214884
f3 = x3 - t^210568
f4 = x4 - t^324307
f5 = x5 - t^320945


I = ideal([f1, f2, f3, f4,f5])

set_verbosity_level(:groebner_walk, 1)

o_t = weight_ordering([1,0,0,0,0,0], degrevlex(R))
o_s = weight_ordering([0,1,1,1,1,1], degrevlex(R))
Ginit = groebner_basis(I, ordering = o_s, complete_reduction = true)
tg = @elapsed Gg = groebner_walk(I, o_t, o_s, algorithm =:generic) #60-70s

ts = @elapsed Gs = groebner_walk(I, o_t, o_s, algorithm =:standard) 

tp = @elapsed Gp = groebner_walk(I, o_t, o_s, algorithm =:perturbed)
tb = @elapsed Gb = groebner_basis(I, ordering = o_t, complete_reduction = true) #140-160s

tf = @elapsed Gf = groebner_basis(I, ordering = o_t, complete_reduction = true, algorithm =:fglm) #error: dimension of ideal must be zero! 

th = @elapsed Gh = groebner_basis(I, ordering = o_t, complete_reduction = true, algorithm =:hilbert)  #throws another error

