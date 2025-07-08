
@doc raw"""
    AffineAlgebraicSet

An affine, geometrically reduced subscheme of an affine space over a field.

# Examples

```jldoctest
julia> A = affine_space(QQ, [:x,:y]);

julia> (x, y) = coordinates(A);

julia> X = algebraic_set(ideal([y - x^2]))
Affine algebraic set
  in 𝔸² over QQ with coordinates x, y
defined by ideal(-x^2 + y)

julia> Y = algebraic_set(ideal([y]))
Affine algebraic set
  in 𝔸² over QQ with coordinates x, y
defined by ideal(y)

julia> Z = set_theoretic_intersection(X, Y)
Affine algebraic set
  in 𝔸² over QQ with coordinates x, y
defined by ideal(-x^2 + y, y)

julia> is_reduced(Z)
true

julia> Z
Affine algebraic set
  in affine 2-space over QQ with coordinates x, y
defined by ideal(y, x)
```
"""

@attributes mutable struct AffineAlgebraicSet{BaseRing<:Field, RingType<:MPolyAnyRing} <: AbsAffineAlgebraicSet{BaseRing, RingType}
  X::AffineScheme
  Xred::AffineScheme
  function AffineAlgebraicSet(X::AffineScheme{S,T}; is_reduced::Bool=false, check::Bool=true) where {S<:Ring,T<:Ring}
    A = new{S, T}()
    A.X = X
    if is_reduced
      A.Xred = X
      set_attribute!(A.Xred, :is_geometrically_reduced, true)
      set_attribute!(A.Xred, :is_reduced, true)
    end
    set_attribute!(A, :is_reduced, true)
    @check is_geometrically_reduced(A) "Algebraic sets must be geometrically reduced"
    set_attribute!(A, :is_geometrically_reduced, true)
    return A
  end
end


