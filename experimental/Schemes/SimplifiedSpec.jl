@attributes mutable struct SimplifiedSpec{BaseRingType, RingType} <: AbsSpec{BaseRingType, RingType} 
  X::AbsSpec
  Y::AbsSpec
  f::AbsSpecMor
  g::AbsSpecMor

  function SimplifiedSpec(X::AbsSpec, Y::AbsSpec, f::AbsSpecMor, g::AbsSpecMor;
      check::Bool=true
    )
    domain(f) === X || error("map is not compatible")
    codomain(f) === Y || error("map is not compatible")
    domain(g) === Y || error("map is not compatible")
    codomain(g) === X || error("map is not compatible")

    if check
      is_identity_map(compose(f, g)) && is_identity_map(compose(g, f)) || error("maps are not inverse to each other")
    end

    result = new{typeof(base_ring(X)), typeof(OO(X))}(X, Y)
    # We need to rewrap the identification maps so that the (co-)domains match
    fwrap = SpecMor(result, Y, pullback(f))
    gwrap = SpecMor(Y, result, pullback(g))
    set_attribute!(fwrap, :inverse, gwrap)
    set_attribute!(gwrap, :inverse, fwrap)
    result.f = fwrap
    result.g = gwrap
    return result 
  end
end

### essential getters
underlying_scheme(X::SimplifiedSpec) = X.X
original(X::SimplifiedSpec) = X.Y
identification_maps(X::SimplifiedSpec) = (X.f, X.g)

### High-level constructors
@Markdown.doc """
    simplify(X::AbsSpec{<:Field})

Given an affine scheme ``X`` with coordinate ring ``R = 𝕜[x₁,…,xₙ]/I`` 
(or a localization thereof), use `Singular`'s `elimpart` to try 
to eliminate variables ``xᵢ`` to arrive at a simpler presentation 
``R ≅ R' = 𝕜[y₁,…,yₘ]/J`` for some ideal ``J``; return 
a `SimplifiedSpec` ``Y`` with ``X`` as its `original`.

***Note:*** The `ambient_coordinate_ring` of the output `Y` will be different
from the one of `X` and hence the two schemes will not compare using `==`.
"""
function simplify(X::AbsSpec{<:Field})
  L, f, g = simplify(OO(X))
  Y = Spec(L)
  YtoX = SpecMor(Y, X, f)
  XtoY = SpecMor(X, Y, g)
  set_attribute!(YtoX, :inverse, XtoY)
  set_attribute!(XtoY, :inverse, YtoX)
  return SimplifiedSpec(Y, X, YtoX, XtoY, check=false)
end

### Methods to roam in the ancestry tree
function some_ancestor(P::Function, X::SimplifiedSpec)
  return P(X) || some_ancestor(P, original(X))
end

function some_ancestor(P::Function, X::PrincipalOpenSubset)
  return P(X) || some_ancestor(P, ambient_scheme(X))
end

@Markdown.doc """
    some_ancestor(P::Function, X::AbsSpec)

Check whether property `P` holds for `X` or some ancestor of `X` in 
case it is a `PrincipalOpenSubset`, or a `SimplifiedSpec`.
"""
function some_ancestor(P::Function, X::AbsSpec)
  return P(X) # This case will only be called when we reached the root.
end

#=
# This crawls up the tree of charts until hitting one of the patches of C.
# Then it returns a pair (f, d) where f is the inclusion morphism 
# of U into the patch V of C and d is a Vector of elements of OO(V)
# which have to be inverted to arrive at U. That is: f induces an 
# isomorphism on the complement of d. 
=#
function _find_chart(U::AbsSpec, C::Covering;
    complement_equations::Vector{T}=elem_type(OO(U))[]
  ) where {T<:RingElem}
  any(W->(W === U), patches(C)) || error("patch not found")
  return identity_map(U), complement_equations
end

function _find_chart(U::PrincipalOpenSubset, C::Covering;
    complement_equations::Vector{T}=elem_type(OO(U))[]
  ) where {T<:RingElem}
  any(W->(W === U), patches(C)) && return identity_map(U), complement_equations
  V = ambient_scheme(U)
  ceq = push!(
              OO(V).(lifted_numerator.(complement_equations)),
              OO(V)(lifted_numerator(complement_equation(U)))
             )
  (f, d) = _find_chart(V, C, complement_equations=ceq)
  return compose(inclusion_morphism(U), f), d
end

function _find_chart(U::SimplifiedSpec, C::Covering;
    complement_equations::Vector{T}=elem_type(OO(U))[]
  ) where {T<:RingElem}
  any(W->(W === U), patches(C)) && return identity_map(U), complement_equations
  V = original(U)
  f, g = identification_maps(U)
  ceq = pullback(g).(complement_equations)
  h, d = _find_chart(V, C, complement_equations=ceq)
  return compose(f, h), d
end

#=
# This follows U in its ancestor tree up to the point 
# where a patch W in C is found. Then it recreates U as a 
# PrincipalOpenSubset UU of W and returns the identification 
# with UU.
=#
function _flatten_open_subscheme(U::AbsSpec, C::Covering)
  any(W->(W === U), patches(C)) || error("patch not found")
  UU = PrincipalOpenSubset(U, one(OO(U)))
  f = inclusion_morphism(UU, U)
  finv = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false))
  set_attribute!(f, :inverse, finv)
  set_attribute!(finv, :inverse, f)
  return f
end

function _flatten_open_subscheme(U::PrincipalOpenSubset, C::Covering)
  some_ancestor(W->any(WW->(WW === W), patches(C)), U) || error("patch not found")
  # TODO: finish implementation
end
