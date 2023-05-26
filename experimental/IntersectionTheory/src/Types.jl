const MPolyDecRingOrQuo = Union{MPolyDecRing, MPolyQuoRing{<:MPolyDecRingElem}}

const MPolyDecRingOrQuoElem = Union{MPolyDecRingElem, MPolyQuoRingElem}

abstract type Bundle end
Base.parent(F::Bundle) = F.parent
rank(F::Bundle) = F.rank

abstract type Variety end
dim(X::Variety) = X.dim

abstract type VarietyHom end
domain(X::VarietyHom) = X.domain
codomain(X::VarietyHom) = X.codomain

Base.show(io::IO, F::Bundle) = print(io, "$(typeof(F).name.name) of rank $(F.rank) on $(F.parent)")
Base.show(io::IO, X::Variety) = print(io, "$(typeof(X).name.name) of dim $(X.dim)")
Base.show(io::IO, f::VarietyHom) = print(io, "$(typeof(f).name.name) from $(f.domain) to $(f.codomain)")

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
    AbstractVarietyMap(X::AbstractVariety, Y::AbstractVariety, fˣ::AffAlgHom, fₓ)
    AbstractVarietyMap(X::AbstractVariety, Y::AbstractVariety, fˣ::Vector, fₓ)
The type of an abstract abstract_variety morphism.
"""
mutable struct AbstractVarietyMap{V1 <: AbstractVarietyT, V2 <: AbstractVarietyT} <: VarietyHom
  domain::V1
  codomain::V2
  dim::Int
  pullback::AffAlgHom
  pushforward::FunctionalMap
  O1::MPolyDecRingOrQuoElem
  T::AbstractBundle{V1}
  function AbstractVarietyMap(X::V1, Y::V2, fˣ::AffAlgHom, fₓ=nothing) where {V1 <: AbstractVarietyT, V2 <: AbstractVarietyT}
    if !(fₓ isa FunctionalMap) && isdefined(X, :point) && isdefined(Y, :point)
      # pushforward can be deduced from pullback in the following cases
      # - explicitly specified (f is relatively algebraic)
      # - X is a point
      # - Y is a point or a curve
      # - all algebraic classes for Y are known
      f_is_alg = fₓ == :alg || dim(X) == 0 || dim(Y) ≤ 1 || get_attribute(Y, :alg) == true
      fₓ = x -> (
	if !f_is_alg
	  @warn "assuming that all algebraic classes are known for\n$Y\notherwise the result may be wrong"
	end;
	sum(integral(xi*fˣ(yi))*di for (i, xi) in zip(dim(Y):-1:0, x[dim(X)-dim(Y):dim(X)])
	    if xi !=0 for (yi, di) in zip(basis(Y, i), dual_basis(Y, i))))
      fₓ = map_from_func(fₓ, X.ring, Y.ring)
    end
    f = new{V1, V2}(X, Y, X.dim-Y.dim, fˣ)
    try
      f.pushforward = fₓ
    catch
    end
    if isdefined(X, :T) && isdefined(Y, :T)
      f.T = AbstractBundle(X, chern_character(X.T) - fˣ(chern_character(Y.T)))
    end
    return f
  end
  function AbstractVarietyMap(X::V1, Y::V2, l::Vector, fₓ=nothing) where {V1 <: AbstractVarietyT, V2 <: AbstractVarietyT}
    # TODO: this fails with check = false
    fˣ = hom(Y.ring, X.ring, l, check = false)
    AbstractVarietyMap(X, Y, fˣ, fₓ)
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
  struct_map::AbstractVarietyMap

  function AbstractVariety(n::Int, R::MPolyDecRingOrQuo)
    base = R isa MPolyQuoRing ? base_ring(base_ring(R)) : base_ring(R)
    X = new(n, R, base)
    set_attribute!(R, :abstract_variety, X)
    set_attribute!(R, :abstract_variety_dim, n)
    return X
  end
end

