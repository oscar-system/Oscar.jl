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

abstract type AbsVarietyT <: Variety end

@doc Markdown.doc"""
    AbsBundle(X::AbsVariety, ch::MPolyDecRingElem)
    AbsBundle(X::AbsVariety, r, c::MPolyDecRingElem)
The type of an abstract bundle.
"""
mutable struct AbsBundle{V <: AbsVarietyT} <: Bundle
  parent::V
  rank::RingElement
  ch::MPolyDecRingOrQuoElem
  chern::MPolyDecRingOrQuoElem
  function AbsBundle(X::V, ch::MPolyDecRingOrQuoElem) where V <: AbsVarietyT
    ch = simplify(ch)
    r = constant_coefficient(ch.f)
    try r = Int(ZZ(QQ(r)))
    catch # r can contain symbolic variables
    end
    new{V}(X, r, ch)
  end
  function AbsBundle(X::V, r::RingElement, c::MPolyDecRingOrQuoElem) where V <: AbsVarietyT
    F = new{V}(X, r)
    F.chern = c
    return F
  end
end

@doc Markdown.doc"""
    AbsVarietyHom(X::AbsVariety, Y::AbsVariety, fˣ::AffAlgHom, fₓ)
    AbsVarietyHom(X::AbsVariety, Y::AbsVariety, fˣ::Vector, fₓ)
The type of an abstract variety morphism.
"""
mutable struct AbsVarietyHom{V1 <: AbsVarietyT, V2 <: AbsVarietyT} <: VarietyHom
  domain::V1
  codomain::V2
  dim::Int
  pullback::AffAlgHom
  pushforward::FunctionalMap
  O1::MPolyDecRingOrQuoElem
  T::AbsBundle{V1}
  function AbsVarietyHom(X::V1, Y::V2, fˣ::AffAlgHom, fₓ=nothing) where {V1 <: AbsVarietyT, V2 <: AbsVarietyT}
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
	    if xi !=0 for (yi, di) in zip(basis(i, Y), dual_basis(i, Y))))
      fₓ = map_from_func(fₓ, X.ring, Y.ring)
    end
    f = new{V1, V2}(X, Y, X.dim-Y.dim, fˣ)
    try
      f.pushforward = fₓ
    catch
    end
    if isdefined(X, :T) && isdefined(Y, :T)
      f.T = AbsBundle(X, ch(X.T) - fˣ(ch(Y.T)))
    end
    return f
  end
  function AbsVarietyHom(X::V1, Y::V2, l::Vector, fₓ=nothing) where {V1 <: AbsVarietyT, V2 <: AbsVarietyT}
    # TODO: this fails with check = false
    fˣ = hom(Y.ring, X.ring, l, check = false)
    AbsVarietyHom(X, Y, fˣ, fₓ)
  end
end


@attributes mutable struct AbsVariety <: AbsVarietyT
  dim::Int
  ring::MPolyDecRingOrQuo
  base::Ring
  point::MPolyDecRingOrQuoElem
  O1::MPolyDecRingOrQuoElem
  T::AbsBundle
  bundles::Vector{AbsBundle}
  struct_map::AbsVarietyHom

  function AbsVariety(n::Int, R::MPolyDecRingOrQuo)
    base = R isa MPolyQuoRing ? base_ring(base_ring(R)) : base_ring(R)
    X = new(n, R, base)
    set_attribute!(R, :variety, X)
    set_attribute!(R, :variety_dim, n)
    return X
  end
end

