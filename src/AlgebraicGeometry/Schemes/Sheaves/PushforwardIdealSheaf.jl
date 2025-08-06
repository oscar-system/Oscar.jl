########################################################################
# PushforwardIdealSheaf                                                #
########################################################################
@attributes mutable struct PushforwardIdealSheaf{SpaceType, OpenType, OutputType,
                                                 RestrictionType
                                                } <: AbsIdealSheaf{
                                                                   SpaceType, OpenType,
                                                                   OutputType, RestrictionType
                                                                  }
  f::CoveredClosedEmbedding
  orig::AbsIdealSheaf
  Ipre::PreSheafOnScheme

  function PushforwardIdealSheaf(
      f::CoveredClosedEmbedding,
      orig::AbsIdealSheaf
    )
    X = domain(f)
    Y = codomain(f)
    @assert X === scheme(orig)

    Ipre = PreSheafOnScheme(Y,
                      OpenType=AbsAffineScheme, OutputType=Ideal,
                      RestrictionType=Map,
                      is_open_func=_is_open_func_for_schemes_without_affine_scheme_open_subscheme(X)
                     )
    I = new{typeof(Y), AbsAffineScheme, Ideal, Map}(f, orig, Ipre)
    return I
  end
end

