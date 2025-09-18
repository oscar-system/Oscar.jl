#=
    Newell's teapot, patch 1
    A 2-dimensional ideal from Tran. "Efficient Groebner walk conversion for implicitization of geometric objects" (2004)
=#
using Oscar

I, o1, o2 = Oscar.newell_patch_with_orderings(QQ, 1)

set_verbosity_level(:groebner_walk, 1)
G = groebner_walk(I, o2, o1; algorithm=:standard)

