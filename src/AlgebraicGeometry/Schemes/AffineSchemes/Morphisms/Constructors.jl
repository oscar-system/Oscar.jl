


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
ERROR: MethodError: no method matching AffineAlgebraicSet(::Spec{QQField, QQMPolyRing}; is_reduced=false, check=false)
Stacktrace:
 [1] algebraic_set(X::Spec{QQField, QQMPolyRing}; is_reduced::Bool, check::Bool)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineAlgebraicSet/Objects/Constructors.jl:15
 [2] variety(X::Spec{QQField, QQMPolyRing}; is_reduced::Bool, check::Bool)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineVariety/Objects/Constructors.jl:12
 [3] affine_space(kk::QQField, n::Int64; variable_name::String)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Constructors.jl:155
 [4] affine_space(kk::QQField, n::Int64)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Constructors.jl:153
 [5] top-level scope
   @ none:1

julia> Y = affine_space(QQ,3)
ERROR: MethodError: no method matching AffineAlgebraicSet(::Spec{QQField, QQMPolyRing}; is_reduced=false, check=false)
Stacktrace:
 [1] algebraic_set(X::Spec{QQField, QQMPolyRing}; is_reduced::Bool, check::Bool)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineAlgebraicSet/Objects/Constructors.jl:15
 [2] variety(X::Spec{QQField, QQMPolyRing}; is_reduced::Bool, check::Bool)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineVariety/Objects/Constructors.jl:12
 [3] affine_space(kk::QQField, n::Int64; variable_name::String)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Constructors.jl:155
 [4] affine_space(kk::QQField, n::Int64)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Constructors.jl:153
 [5] top-level scope
   @ none:1

julia> SpecMor(X, Y, gens(OO(X)));
ERROR: UndefVarError: X not defined
Stacktrace:
 [1] top-level scope
   @ none:1
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
  with coordinates x1 x2 x3
  over rational field

julia> identity_map(X);
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
  with coordinates x1 x2 x3
  over rational field

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of multivariate polynomial ring by ideal with 1 generator

julia> f = inclusion_morphism(Y, X);

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
  with coordinates x1 x2 x3
  over rational field

julia> R = OO(X)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> Y = subscheme(X, x1)
Spec of Quotient of multivariate polynomial ring by ideal with 1 generator

julia> m1 = inclusion_morphism(Y, X);

julia> m2 = identity_map(X);

julia> compose(m1, m2) == m1
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
ERROR: MethodError: no method matching AffineAlgebraicSet(::Spec{QQField, QQMPolyRing}; is_reduced=false, check=false)
Stacktrace:
 [1] algebraic_set(X::Spec{QQField, QQMPolyRing}; is_reduced::Bool, check::Bool)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineAlgebraicSet/Objects/Constructors.jl:15
 [2] variety(X::Spec{QQField, QQMPolyRing}; is_reduced::Bool, check::Bool)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineVariety/Objects/Constructors.jl:12
 [3] affine_space(kk::QQField, n::Int64; variable_name::String)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Constructors.jl:155
 [4] affine_space(kk::QQField, n::Int64)
   @ Oscar ~/Kaiserslautern/neuer_Oscar_Klon/Oscar.jl/src/AlgebraicGeometry/Schemes/AffineSchemes/Objects/Constructors.jl:153
 [5] top-level scope
   @ none:1

julia> R = OO(X)
ERROR: UndefVarError: X not defined
Stacktrace:
 [1] top-level scope
   @ none:1

julia> (x1,x2,x3) = gens(R)
ERROR: UndefVarError: R not defined
Stacktrace:
 [1] top-level scope
   @ none:1

julia> Y = subscheme(X, x1)
ERROR: UndefVarError: X not defined
Stacktrace:
 [1] top-level scope
   @ none:1

julia> restrict(identity_map(X), Y, Y) == identity_map(Y)
ERROR: UndefVarError: X not defined
Stacktrace:
 [1] top-level scope
   @ none:1
```
"""
function restrict(f::SpecMor, U::AbsSpec, V::AbsSpec; check::Bool=true)
  @check issubset(U, domain(f)) "second argument does not lie in the domain of the map"
  @check issubset(V, codomain(f)) "third argument does not lie in the codomain of the map"
  @check issubset(U, preimage(f, V)) "the image of the restriction is not contained in the restricted codomain"
  imgs = pullback(f).(gens(domain(pullback(f))))
  return SpecMor(U, V, OO(U).(imgs), check=check)
end
