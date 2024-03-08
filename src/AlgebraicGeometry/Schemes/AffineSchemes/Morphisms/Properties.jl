

########################################
# (1) Isomorphism, inverse and identity
########################################

@doc raw"""
    is_isomorphism(f::AbsAffineSchemeMor)

This method checks if a morphism is an isomorphism.
"""
@attr Bool function is_isomorphism(f::AbsAffineSchemeMor)
  has_attribute(f, :inverse) && return true
  is_isomorphism(pullback(f)) || return false
  set_attribute!(f, :inverse, morphism(codomain(f), domain(f), inverse(pullback(f))))
  return true
end

@doc raw"""
    is_inverse_of(f::AbsAffineSchemeMor, g::AbsAffineSchemeMor)

This method checks if a morphism ``f`` is the inverse of a morphism ``g``.
"""
function is_inverse_of(f::S, g::T) where {S<:AbsAffineSchemeMor, T<:AbsAffineSchemeMor}
  return is_isomorphism(f) && (inverse(f) == g)
end

@doc raw"""
    is_identity_map(f::AbsAffineSchemeMor)

This method checks if a morphism is the identity map.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> is_identity_map(inclusion_morphism(Y, X))
false
```
"""
is_identity_map(f::AbsAffineSchemeMor) = (domain(f) == codomain(f)) && all(x->(pullback(f)(x) == x), gens(OO(domain(f))))

@attr Bool function is_birational(f::AbsAffineSchemeMor)
  error("verification of birationality not implemented")
end

