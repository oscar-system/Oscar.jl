

########################################
# (1) Isomorphism, inverse and identity
########################################

@doc raw"""
    is_isomorphism(f::AbsSpecMor)

This method checks if a morphism is an isomorphism.
"""
@attr Bool function is_isomorphism(f::AbsSpecMor)
  has_attribute(f, :inverse) && return true
  is_isomorphism(pullback(f)) || return false
  set_attribute!(f, :inverse, SpecMor(codomain(f), domain(f), inverse(pullback(f))))
  return true
end

@doc raw"""
    is_inverse_of(f::AbsSpecMor, g::AbsSpecMor)

This method checks if a morphism ``f`` is the inverse of a morphism ``g``.
"""
function is_inverse_of(f::S, g::T) where {S<:AbsSpecMor, T<:AbsSpecMor}
  return is_isomorphism(f) && (inverse(f) == g)
end

@doc raw"""
    is_identity_map(f::AbsSpecMor)

This method checks if a morphism is the identity map.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine 3-space
 over Rational Field
with coordinates
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> is_identity_map(inclusion_morphism(Y, X))
false
```
"""
is_identity_map(f::AbsSpecMor) = (domain(f) == codomain(f)) && all(x->(pullback(f)(x) == x), gens(OO(domain(f))))
