################################################################################
#
#  Tropicalization of affine linear ideals
#
#  WARNING: assumes without test that `gens(I)` has degree one
#
################################################################################

function tropical_variety_affine_linear(I::MPolyIdeal,nu::TropicalSemiringMap; weighted_polyhedral_complex_only::Bool=false)
    # compute a reduced GB to check whether I is homogeneous
    G = groebner_basis(I,complete_reduction=true)

    if all(Oscar._is_homogeneous,G)
        # input homogneeous, construct TropicalVariety via TropicalLinearSpace
        return tropical_variety(tropical_linear_space(I,nu,weighted_polyhedral_complex_only=weighted_polyhedral_complex_only))
    end

    # input inhomogeneous, homogenise first
    Ih,_,_ = homogenize_pre_tropicalization(I)
    TropLh = tropical_linear_space(Ih,nu,weighted_polyhedral_complex_only=true)
    return dehomogenize_post_tropicalization(TropLh)
end
