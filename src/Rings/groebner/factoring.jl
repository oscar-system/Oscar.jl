@doc raw"""
    factoring_standard_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I))) 

Return a factorization of the standard basis of `I` with respect to `ordering`,
i.e. a list of IdealGens.
The intersection of these ideals has the same zero set as the input ideal, i.e.
the radical of the intersection coincides with the radical of the input ideal.
"""
function factoring_standard_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)))

    #= apply highest corner standard basis variant in Singular =#
    sL = Singular.facstd(singular_generators(I, ordering))

    L = [Oscar.IdealGens(base_ring(I), r, false) for r in sL]

    return L
end
