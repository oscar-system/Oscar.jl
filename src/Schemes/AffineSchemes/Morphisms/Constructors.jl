export SpecMor, identity_map, inclusion_morphism, restrict, compose



########################################################################
# (1) General constructors
########################################################################

@Markdown.doc """
    SpecMor(X::AbsSpec, Y::AbsSpec, f::Vector{<:RingElem}; check::Bool=true)

This method constructs a morphism from the scheme ``X``
to the scheme ``Y``. For this one has to specify the images
of the coordinates (the generators of `ambient_ring(Y)`)
under the pullback map ``𝒪(Y) → 𝒪(X)`` as third argument.

Note that expensive checks can be turned off by setting `check=false`.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> Y = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> SpecMor(X, Y, gens(OO(X)));
```
"""
function SpecMor(
      X::AbsSpec,
      Y::AbsSpec,
      f::Vector{<:RingElem};
      check::Bool=true
  )
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f), check=check), check=check)
end

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec{<:Ring, <:MPolyRing},
      f::Vector{<:RingElem};
      check::Bool=true
  )
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f)), check=check)
end

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec,
      f::Vector;
      check::Bool=true
  )
  return SpecMor(X, Y, OO(X).(f), check=check)
end



########################################################################
# (2) Special constructors
########################################################################

@Markdown.doc """
    identity_map(X::AbsSpec{<:Any, <:MPolyRing})

This method constructs the identity morphism from an affine scheme to itself.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> identity_map(X);
```
"""
identity_map(X::AbsSpec{<:Any, <:MPolyRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(OO(X))))
identity_map(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X)), check=false))
identity_map(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X)), check=false))
identity_map(X::AbsSpec{<:Any, <:MPolyQuo}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X))))


@Markdown.doc """
    inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true)

This method constructs the inclusion map from ``X`` to ``Y``.
For convenience, also the method `inclusion_morphism` is supported.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> inclusion_morphism(Y, X);
```
"""
inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true) = SpecMor(X, Y, gens(ambient_ring(Y)), check=check)


@Markdown.doc """
    compose(f::AbsSpecMor, g::AbsSpecMor)

This method computes the composition of two morphisms.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> m1 = inclusion_morphism(Y, X);

julia> m2 = identity_map(X);

julia> compose(m1, m2) == m1
true
```
"""
function compose(f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)
  codomain(f) == domain(g) || error("Morphisms can not be composed")
  return SpecMor(domain(f), codomain(g), compose(pullback(g), pullback(f)), check=check)
end


@Markdown.doc """
    restrict(f::SpecMor, U::AbsSpec, V::AbsSpec)

This method restricts the domain of the morphism ``f``
to ``U`` and its codomain to ``V``.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Spec of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> R = OO(X)
Multivariate Polynomial Ring in x1, x2, x3 over Rational Field

julia> (x1,x2,x3) = gens(R)
3-element Vector{fmpq_mpoly}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of Multivariate Polynomial Ring in x1, x2, x3 over Rational Field by ideal(x1)

julia> restrict(identity_map(X), Y, Y) == identity_map(Y)
true
```
"""
function restrict(f::SpecMor, U::AbsSpec, V::AbsSpec; check::Bool=true)
  if check
    issubset(U, domain(f)) || error("second argument does not lie in the domain of the map")
    issubset(V, codomain(f)) || error("third argument does not lie in the codomain of the map")
    issubset(U, preimage(f, V)) || error("the image of the restriction is not contained in the restricted codomain")
  end
  return SpecMor(U, V, OO(U).(pullback(f).(gens(domain(pullback(f))))), check=check)
end
