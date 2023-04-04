########################################################
# (1) Display
########################################################

function Base.show(io::IO, X::ProjectiveAlgebraicSet)
  d = dim(ambient_space(X))
  println(io, "Projective algebraic set in $(ambient_space(X))")
end

########################################################
# (2) Irreducible Components
########################################################

@doc Markdown.doc"""
    irreducible_components(X::ProjectiveAlgebraicSet) -> Vector{ProjectiveVariety}

Return the irreducible components of `X` defined over the same base field.

Note that even if `X` is irreducible, there may be several geometric irreducible components.
"""
function irreducible_components(X::ProjectiveAlgebraicSet)
  throw(NotImplementedError())
end


@doc Markdown.doc"""
    geometric_irreducible_components(X::ProjectiveAlgebraicSet) -> Vector{ProjectiveVariety}

Return the geometric irreducible components of `X`.

They are the irreducible components of `X` seen over an algebraically closed field.

This is expensive and involves taking field extensions.
"""
function geometric_irreducible_components(X::ProjectiveAlgebraicSet)
  throw(NotImplementedError())
end
