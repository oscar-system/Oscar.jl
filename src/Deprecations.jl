# Deprecated in 0.10.*

@deprecate automorphisms(x::Graph) automorphism_group_generators(x)
@deprecate _containement_helper(R::MPolyRing, N::Int, M::Int, I::MPolyIdeal, W::Vector, ord::Symbol) _containment_helper(R, N, M, I, W, ord)
