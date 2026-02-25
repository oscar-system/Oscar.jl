const MPolyDecRingOrQuo = Union{MPolyDecRing, MPolyQuoRing{<:MPolyDecRingElem}}

const MPolyDecRingOrQuoElem = Union{MPolyDecRingElem, MPolyQuoRingElem}

# TODO Bundle or AbstractBundle?
abstract type Bundle end


@doc raw"""
    parent(F::AbstractBundle)

Return the underlying abstract variety of `F`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 6

julia> parent(Q) == G
true

```
"""
Base.parent(F::Bundle) = F.parent

@doc raw"""
    rank(F::AbstractBundle)

Return the rank of `F`.

# Examples
```jldoctest
julia> G = abstract_grassmannian(3, 5)
AbstractVariety of dim 6

julia> Q = tautological_bundles(G)[2]
AbstractBundle of rank 2 on AbstractVariety of dim 6

julia> rank(symmetric_power(Q, 3))
4

```
"""
rank(F::Bundle) = F.rank

abstract type Variety end

@doc raw"""
    dim(X::AbstractVariety)

Return the dimension of `X`.

# Examples
```jldoctest
julia> P2 = abstract_projective_space(2)
AbstractVariety of dim 2

julia> P3 = abstract_projective_space(3)
AbstractVariety of dim 3

julia> dim(P2*P3)
5

```
"""
dim(X::Variety) = X.dim

abstract type VarietyHom end

@doc raw"""
    domain(f::AbstractVarietyMap)

Return the domain of `f`.
"""
domain(X::VarietyHom) = X.domain

@doc raw"""
    codomain(f::AbstractVarietyMap)

Return the codomain of `f`.
"""
codomain(X::VarietyHom) = X.codomain

function Base.show(io::IO, F::Bundle)
  if Oscar.is_terse(io)
    print(io, "AbstractBundle")
  else
    print(io, "AbstractBundle of rank $(F.rank) on $(F.parent)")
  end
end

function Base.show(io::IO, X::Variety)
  if Oscar.is_terse(io)
    print(io, "AbstractVariety")
  else
    print(io, "AbstractVariety of dim $(X.dim)")
  end
end

function Base.show(io::IO, f::VarietyHom)
  if Oscar.is_terse(io)
    print(io, "AbstractVarietyMap")
  else
    print(io, "AbstractVarietyMap from $(f.domain) to $(f.codomain)")
  end
end

abstract type AbstractVarietyT <: Variety end

@doc raw"""
    AbstractBundle(X::AbstractVariety, ch::MPolyDecRingElem)
    AbstractBundle(X::AbstractVariety, r, c::MPolyDecRingElem)

The type of an abstract bundle.
"""
mutable struct AbstractBundle{V <: AbstractVarietyT} <: Bundle
  parent::V
  rank::RingElement
  ch::MPolyDecRingOrQuoElem
  chern::MPolyDecRingOrQuoElem

  function AbstractBundle(X::V, ch::MPolyDecRingOrQuoElem) where V <: AbstractVarietyT
    ch = simplify(ch)
    r = constant_coefficient(ch.f)
    try r = Int(ZZ(QQ(r)))
    catch # r can contain symbolic variables
    end
    new{V}(X, r, ch)
  end

  function AbstractBundle(X::V, r::RingElement, c::MPolyDecRingOrQuoElem) where V <: AbstractVarietyT
    F = new{V}(X, r)
    F.chern = c
    return F
  end
end

@doc raw"""
    AbstractVarietyMap(X::AbstractVariety, Y::AbstractVariety, f_pullback::AffAlgHom, f_pushforward)
    AbstractVarietyMap(X::AbstractVariety, Y::AbstractVariety, f_pullback::Vector, f_pushforward)

The type of an abstract variety morphism.
"""
mutable struct AbstractVarietyMap{V1 <: AbstractVarietyT, V2 <: AbstractVarietyT} <: VarietyHom
  domain::V1
  codomain::V2
  dim::Int
  pullback::AffAlgHom
  pushforward::MapFromFunc
  O1::MPolyDecRingOrQuoElem
  T::AbstractBundle{V1}

  function AbstractVarietyMap(X::V1, Y::V2, f_pullback::AffAlgHom, f_pushforward=nothing) where {V1 <: AbstractVarietyT, V2 <: AbstractVarietyT}
    if !(f_pushforward isa MapFromFunc) && isdefined(X, :point) && isdefined(Y, :point)
      # pushforward can be deduced from pullback in the following cases
      # - explicitly specified (f is relatively algebraic)
      # - X is a point
      # - Y is a point or a curve
      # - all algebraic classes for Y are known
      f_is_alg = f_pushforward == :alg || dim(X) == 0 || dim(Y) ≤ 1 || get_attribute(Y, :alg) == true
      f_pushforward = x -> (
        if !f_is_alg
          @warn "assuming that all algebraic classes are known for\n$Y\notherwise the result may be wrong"
        end;
        sum(integral(xi*f_pullback(yi))*di for (i, xi) in zip(dim(Y):-1:0, x[dim(X)-dim(Y):dim(X)])
            if xi !=0 for (yi, di) in zip(basis(Y, i), dual_basis(Y, i))))
      f_pushforward = MapFromFunc(X.ring, Y.ring, f_pushforward)
    end
    f = new{V1, V2}(X, Y, X.dim-Y.dim, f_pullback)
    try
      f.pushforward = f_pushforward
    catch
    end
    if isdefined(X, :T) && isdefined(Y, :T)
      f.T = AbstractBundle(X, chern_character(X.T) - f_pullback(chern_character(Y.T)))
    end
    return f
  end

  function AbstractVarietyMap(X::V1, Y::V2, l::Vector, f_pushforward=nothing) where {V1 <: AbstractVarietyT, V2 <: AbstractVarietyT}
    # TODO: this fails with check = false
    f_pullback = hom(Y.ring, X.ring, l, check = false)
    AbstractVarietyMap(X, Y, f_pullback, f_pushforward)
  end
end


@attributes mutable struct AbstractVariety <: AbstractVarietyT
  dim::Int
  ring::MPolyDecRingOrQuo
  base::Ring
  point::MPolyDecRingOrQuoElem
  O1::MPolyDecRingOrQuoElem
  T::AbstractBundle
  bundles::Vector{AbstractBundle}
  structure_map::AbstractVarietyMap

  function AbstractVariety(n::Int, R::MPolyDecRingOrQuo)
    base = R isa MPolyQuoRing ? base_ring(base_ring(R)) : base_ring(R)
    X = new(n, R, base)
    set_attribute!(R, :abstract_variety, X)
    set_attribute!(R, :abstract_variety_dim, n)
    return X
  end
end
