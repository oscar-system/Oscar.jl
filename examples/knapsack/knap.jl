using Oscar 
using GroebnerWalk
R, (t, x1, x2, x3, x4, x5, x6) = polynomial_ring(QQ, ["t","x1", "x2", "x3", "x4", "x5", "x6"])






f1 = x1 - t^53725
f2 = x2 - t^61919
f3 = x3 - t^64670
f4 = x4 - t^69340
f5 = x5 - t^78539
f6 = x6 - t^95043 

I = ideal([f1, f2, f3, f4, f5, f6])

set_verbosity_level(:groebner_walk, 1)

o_t = weight_ordering([1,0,0,0,0,0,0], degrevlex(R))
o_s = weight_ordering([0,1,1,1,1,1,1], degrevlex(R))
Ginit = groebner_basis(I, ordering = o_s, complete_reduction = true)
tg = @elapsed Gg = groebner_walk(I, o_t, o_s, algorithm =:generic) #300s

ts = @elapsed Gs = groebner_walk(I, o_t, o_s, algorithm =:standard) #throws an error...oops

tp = @elapsed Gp = groebner_walk(I, o_t, o_s, algorithm =:perturbed) #also throws an error
tb = @elapsed Gb = groebner_basis(I, ordering = o_t, complete_reduction = true) 

tf = @elapsed Gf = groebner_basis(I, ordering = o_t, complete_reduction = true, algorithm =:fglm) #error: dimension of ideal must be zero! 

th = @elapsed Gh = groebner_basis(I, ordering = o_t, complete_reduction = true, algorithm =:hilbert) 