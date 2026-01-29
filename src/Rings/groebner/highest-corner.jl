@doc raw"""
    standard_basis_highest_corner(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I))) 

Return a standard basis of `I` with respect to `ordering`. `ordering` needs to be local, the coefficient ring needs to be `QQ`.
The algorithm first computes a standard basis over a finite field in order to get an upper bound for the highest corner fast.
Then this bound is used to speed up the standard basis computation over `QQ´.
"""
function standard_basis_highest_corner(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)))
    @req is_local(ordering) "Monomial ordering must be local for this variant."
    @req coefficient_ring(I) == QQ "Base ring must be QQ."

    #= apply highest corner standard basis variant in Singular =#
    ssb = Singular.LibStandard.groebner(singular_groebner_generators(I, ordering), "HC")

    sb = IdealGens(base_ring(I), ssb, false)
    sb.isGB = true
    sb.ord  = ordering
    if isdefined(sb, :S)
        sb.S.isGB = true
    end
    I.gb[ordering] = sb
    return sb
end
