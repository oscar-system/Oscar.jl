########################################################################
#
# (1) AbsAffineScheme interface
#
########################################################################

@doc raw"""
    underlying_scheme(X::AffineAlgebraicSet) -> AbsAffineScheme

Return the underlying reduced scheme defining ``X``.

This is used to forward the `AbsAffineScheme` functionality to ``X``, but may
trigger the computation of a radical ideal. Hence this can be expensive.
"""
function underlying_scheme(X::AffineAlgebraicSet)
  if isdefined(X, :Xred)
    return X.Xred
  end
  X.Xred = reduced_scheme(X.X)[1]
  return X.Xred
end

function underlying_scheme(X::AffineAlgebraicSet{<:Field,<:MPolyQuoRing})
  if isdefined(X, :Xred)
    return X.Xred
  end
  # avoid constructing a morphism
  F = fat_scheme(X)
  I = saturated_ideal(defining_ideal(F))
  if has_attribute(F, :is_reduced) && is_reduced(F)
      Irad = I
  else
    Irad = radical(I)
  end
  X.Xred = spec(base_ring(Irad), Irad)
  return X.Xred
end

########################################################################
#
# (2) AbsAffineAlgebraicSet interface
#
########################################################################
@doc raw"""
    fat_scheme(X::AffineAlgebraicSet) -> AbsAffineScheme

Return a scheme whose reduced subscheme is ``X``.

This does not trigger any computation and is therefore cheap.
Use this instead of `underlying_scheme` when possible.
"""
function fat_scheme(X::AffineAlgebraicSet)
  return X.X
end

########################################################################
#
# (3) Further attributes
#
########################################################################

@doc raw"""
    vanishing_ideal(X::AbsAffineAlgebraicSet) -> Ideal

Return the ideal of all polynomials vanishing in ``X``.

By Hilbert's Nullstellensatz this is a radical ideal.
!!! note
    This triggers the computation of a radical, which is expensive.
"""
vanishing_ideal(X::AbsAffineAlgebraicSet) = saturated_ideal(defining_ideal(X))
#vanishing_ideal(X::AbsAffineAlgebraicSet) = ambient_closure_ideal(X) # TODO: What about this here? Should it also be removed?

@doc raw"""
    fat_ideal(X::AbsAffineAlgebraicSet) -> Ideal

Return an ideal whose radical is the vanishing ideal of `X`.

If `X` is constructed from an ideal `I` this returns `I`.

```jldoctest
julia> A2 = affine_space(QQ, [:x,:y])
Affine space of dimension 2
  over rational field
with coordinates [x, y]

julia> (x, y) = coordinates(A2);

julia> I = ideal([x^2, y]);

julia> X = algebraic_set(I)
Affine algebraic set
  in affine 2-space over QQ with coordinates [x, y]
defined by ideal (x^2, y)

julia> fat_ideal(X) === I
true
```
"""
fat_ideal(X::AbsAffineAlgebraicSet) = saturated_ideal(defining_ideal(fat_scheme(X)))

# avoid computing the underlying scheme
ambient_space(X::AbsAffineAlgebraicSet) = ambient_space(fat_scheme(X))

ambient_space(X::AbsAffineAlgebraicSet{S,T}) where {S<:Field, T<:MPolyRing} = X


