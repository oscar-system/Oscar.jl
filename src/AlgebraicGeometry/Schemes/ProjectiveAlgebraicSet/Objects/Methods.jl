########################################################
# (1) Display
########################################################


function Base.show(io::IO, ::MIME"text/plain", X::AbsProjectiveAlgebraicSet{<:Field,<:MPolyQuoRing})
    println(io, "Vanishing locus ")
    println(io, "  in ", ambient_space(X))
    print(io, "  of ")
    Base.show(io, MIME("text/plain"), defining_ideal(X))
end


function Base.show(io::IO, X::AbsProjectiveAlgebraicSet)
    if get(io, :supercompact, false)
      print(io, "Projective algebraic set")
    else
      print(io, "Projective algebraic set in ")
      print(io, ambient_space(X))
    end
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
  I = defining_ideal(X)
  J = minimal_primes(I)
  return typeof(X)[vanishing_locus(j,check=false) for j in J]
end


@doc Markdown.doc"""
    geometric_irreducible_components(X::ProjectiveAlgebraicSet) -> Vector{ProjectiveVariety}

Return the geometric irreducible components of `X`.

They are the irreducible components of `X` seen over an algebraically closed field.

This is expensive and involves taking field extensions.
"""
function geometric_irreducible_components(X::ProjectiveAlgebraicSet)
  throw(NotImplementedError(:ProjectiveAlgebraicSet))
end
