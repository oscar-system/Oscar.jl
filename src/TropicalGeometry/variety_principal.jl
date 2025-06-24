################################################################################
#
#  Tropicalization of principal ideals
#
#  WARNING: assumes without test that `gens(I)` has cardinality one
#
################################################################################

function tropical_variety_principal(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    ###
    # Construct TropicalVariety from TropicalHypersurface
    ###
    g = gen(I,1)
    TropV = tropical_variety(tropical_hypersurface(g,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only))
    if !weighted_polyhedral_complex_only
        set_attribute!(TropV,:algebraic_ideal,I)
        set_attribute!(TropV,:tropical_semiring_map,nu)
    end
    return TropV
end
