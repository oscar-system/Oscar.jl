########################################################
# (1) Display
########################################################

function Base.show(io::IO, X::AffineAlgebraicSet)
  d = dim(ambient_space(X))
  println(io, "affine algebraic set in $(ambient_space(X))")
end

########################################################
# (2) Irreducible Components
########################################################

@doc raw"""
    irreducible_components(X::AffineAlgebraicSet) -> Vector{AffineVariety}

Return the irreducible components of `X` defined over the same base field.

Note that this is not the same as the geometric irreducible components.
"""
function irreducible_components(X::AffineAlgebraicSet)
  throw(NotImplementedError())
end


@doc raw"""
    geometric_irreducible_components(X::AffineAlgebraicSet) -> Vector{AffineVariety}

Return the geometric irreducible components of `X`.

They are the irreducible components of `X` seen over an algebraically closed field.

This is expensive and involves taking field extensions.
"""
function geometric_irreducible_components(X::AffineAlgebraicSet)
  throw(NotImplementedError())
end
