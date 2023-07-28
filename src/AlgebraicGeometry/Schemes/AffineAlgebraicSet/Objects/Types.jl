
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
  X::Spec
  Xred::Spec
  function AffineAlgebraicSet(X::Spec; is_reduced::Bool=false, check::Bool=true)
    A = new{typeof(base_ring(X)), typeof(OO(X))}()
    A.X = X
    if check && is_reduced
      is_geometrically_reduced(X) || error("Algebraic sets must be geometrically reduced")
    end
    if is_reduced
      A.Xred = X
      set_attribute!(A.Xred, :is_geometrically_reduced, true)
      set_attribute!(A.Xred, :is_reduced, true)
      # unlock the scheme methods for geometrically reduced schemes.
    end
    return A
  end
end


