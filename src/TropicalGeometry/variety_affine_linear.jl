################################################################################
#
#  Tropicalization of affine linear ideals
#
#  WARNING: assumes without test that `gens(I)` has degree one
#
################################################################################

function tropical_variety_affine_linear(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    ###
    # Compute reduced Groebner basis (usually already cached),
    # and check whether the linear polynomials have a constant term
    ###
    R = base_ring(I)
    G = groebner_basis(I,complete_reduction=true)
    if min(total_degree.(Iterators.flatten(collect.(terms.(G))))...)==1
        # input homogneeous, construct TropicalVariety via TropicalLinearSpace
        TropV = tropical_variety(tropical_linear_space(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only))
        if !weighted_polyhedral_complex_only
            set_attribute!(TropV,:algebraic_ideal,I)
            set_attribute!(TropV,:tropical_semiring_map,nu)
        end
        return TropV
    else
        # input inhomogeneous, homogenise first
        Ih = homogenize_pre_tropicalization(I)
        TropLh = tropical_linear_space(Ih,nu,weighted_polyhedral_complex_only=true)
        Sigma = dehomogenize_post_tropicalization(polyhedral_complex(TropLh))

        multiplicities = ones(ZZRingElem, n_maximal_polyhedra(Sigma))
        TropV = tropical_variety(Sigma,multiplicities)
        if !weighted_polyhedral_complex_only
            set_attribute!(TropV,:algebraic_ideal,I)
            set_attribute!(TropV,:tropical_semiring_map,nu)
        end
        return TropV
    end
end
