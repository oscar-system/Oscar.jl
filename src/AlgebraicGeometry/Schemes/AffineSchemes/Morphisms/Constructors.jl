


########################################################################
# (1) General constructors
########################################################################

@doc raw"""
    SpecMor(X::AbsSpec, Y::AbsSpec, f::Vector{<:RingElem}; check::Bool=true)

This method constructs a morphism from the scheme ``X``
to the scheme ``Y``. For this one has to specify the images
of the coordinates (the generators of `ambient_coordinate_ring(Y)`)
under the pullback map ``𝒪(Y) → 𝒪(X)`` as third argument.

Note that expensive checks can be turned off by setting `check=false`.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> SpecMor(X, Y, gens(OO(X)))
Affine scheme morphism
  from [x1, x2, x3]  affine 3-space over QQ
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3
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
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f), check=check), check=check)
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

@doc raw"""
    identity_map(X::AbsSpec{<:Any, <:MPolyRing})

This method constructs the identity morphism from an affine scheme to itself.

# Examples
```jldoctest
julia> X = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> identity_map(X)
Affine scheme morphism
  from [x1, x2, x3]  affine 3-space over QQ
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3
```
"""
identity_map(X::AbsSpec{<:Any, <:MPolyRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(OO(X)), check=false), check=false)
identity_map(X::AbsSpec{<:Any, <:MPolyQuoLocRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_coordinate_ring(X)), check=false), check=false)
identity_map(X::AbsSpec{<:Any, <:MPolyLocRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_coordinate_ring(X)), check=false), check=false)
identity_map(X::AbsSpec{<:Any, <:MPolyQuoRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_coordinate_ring(X)), check=false), check=false)


@doc raw"""
    inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true)

Return the inclusion map from ``X`` to ``Y``.

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
    by ideal(x1)

julia> f = inclusion_morphism(Y, X)
Affine scheme morphism
  from [x1, x2, x3]  V(x1)
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> 0
  x2 -> x2
  x3 -> x3

julia> I = kernel(pullback(f))  # this is a way to obtain the ideal ``I ⊆  O(X)`` cutting out ``Y`` from ``X``.
ideal(x1)

julia> base_ring(I) == OO(X)
true
```
"""
inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true) = SpecMor(X, Y, gens(ambient_coordinate_ring(Y)), check=check)


@doc raw"""
    compose(f::AbsSpecMor, g::AbsSpecMor)

This method computes the composition of two morphisms.

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
    by ideal(x1)

julia> m1 = inclusion_morphism(Y, X)
Affine scheme morphism
  from [x1, x2, x3]  V(x1)
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> 0
  x2 -> x2
  x3 -> x3

julia> m2 = identity_map(X)
Affine scheme morphism
  from [x1, x2, x3]  affine 3-space over QQ
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

julia> m3 = identity_map(Y)
Affine scheme morphism
  from [x1, x2, x3]  V(x1)
  to   [x1, x2, x3]  V(x1)
given by the pullback function
  x1 -> 0
  x2 -> x2
  x3 -> x3

julia> compose(m3, compose(m1, m2)) == m1
true
```
"""
function compose(f::AbsSpecMor, g::AbsSpecMor)
  codomain(f) === domain(g) || error("Morphisms can not be composed")
  return SpecMor(domain(f), codomain(g), compose(pullback(g), pullback(f)), check=false)
end


@doc raw"""
    restrict(f::SpecMor, U::AbsSpec, V::AbsSpec)

This method restricts the domain of the morphism ``f``
to ``U`` and its codomain to ``V``.

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
    by ideal(x1)

julia> restrict(identity_map(X), Y, Y) == identity_map(Y)
true
```
"""
function restrict(f::SpecMor, U::AbsSpec, V::AbsSpec; check::Bool=true)
  @check issubset(U, domain(f)) "second argument does not lie in the domain of the map"
  @check issubset(V, codomain(f)) "third argument does not lie in the codomain of the map"
  @check issubset(U, preimage(f, V)) "the image of the restriction is not contained in the restricted codomain"
  return SpecMor(U, V, (x->OO(U)(x, check=check)).(pullback(f).(gens(domain(pullback(f))))), check=check)
end
