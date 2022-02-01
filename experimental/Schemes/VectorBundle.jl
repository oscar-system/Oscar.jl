@attribute mutable struct VectorBundle{
    CoveredSchemeType<:CoveredScheme,
    CoveringType<:Covering,
    SpecType<:Spec,
    MatrixType
  }
  X::CoveredSchemeType
  C::CoveringType
  transitions::Dict{Tuple{SpecType, SpecType}, MatrixType}

  function VectorBundle(
      X::CoveredScheme, 
      C::Covering, 
      g::Dict{Tuple{SpecType, SpecType}, MatrixType}; 
      check::Bool=true
    ) where {SpecType<:Spec, MatrixType}
    C in coverings(X) || error("covering is not listed")
    for (U, V) in keys(g)
      (U in patches(C) && V in patches(C)) || error("affine patches not listed")
    end
    if check
      #TODO: Do some more sophisticated checks
    end
    return new{typeof(X), typeof(C), affine_patch_type(X), MatrixType}(X, C, g)
  end
end

scheme(E::VectorBundle) = E.X
covering(E::VectorBundle) = E.C
transitions(E::VectorBundle) = E.transitions
getindex(E::VectorBundle, U::SpecType, V::SpecType) where {SpecType<:Spec} = E.transitions[(U, V)]

set_name!(E::VectorBundle, name::String) = set_attribute!(E, :name, name)
has_name(E::VectorBundle) = has_attribute(E, :name)
name_of(E::VectorBundle) = get_attribute!(E, :name)::String


