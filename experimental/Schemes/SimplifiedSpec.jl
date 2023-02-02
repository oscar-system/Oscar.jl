########################################################################
# Simplified Spectra                                                   #
########################################################################

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
function has_ancestor(P::Function, X::SimplifiedSpec)
  return P(X) || has_ancestor(P, original(X))
end

function has_ancestor(P::Function, X::PrincipalOpenSubset)
  return P(X) || has_ancestor(P, ambient_scheme(X))
end

@Markdown.doc """
    has_ancestor(P::Function, X::AbsSpec)

Check whether property `P` holds for `X` or some ancestor of `X` in 
case it is a `PrincipalOpenSubset`, or a `SimplifiedSpec`.
"""
function has_ancestor(P::Function, X::AbsSpec)
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
  ceq = Vector{elem_type(OO(V))}(pullback(g).(complement_equations))
  h, d = _find_chart(V, C, complement_equations=ceq)
  return compose(f, h), d
end

function __find_chart(U::AbsSpec, C::Covering) 
  any(W->(W === U), patches(C)) || error("patch not found")
  return U
end

function __find_chart(U::PrincipalOpenSubset, C::Covering)
  any(W->(W === U), patches(C)) && return U
  return __find_chart(ambient_scheme(U), C)
end

function __find_chart(U::SimplifiedSpec, C::Covering)
  any(W->(W === U), patches(C)) && return U
  return __find_chart(original(U), C)
end

#=
# This follows U in its ancestor tree up to the point 
# where a patch W in C is found. Then it recreates U as a 
# PrincipalOpenSubset UU of W and returns the identification 
# with UU.
=#
function _flatten_open_subscheme(
    U::PrincipalOpenSubset, C::Covering;
    iso::AbsSpecMor=begin
      UU = PrincipalOpenSubset(U, one(OO(U)))
      f = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false), check=false)
      f_inv = SpecMor(UU, U, hom(OO(U), OO(UU), gens(OO(UU)), check=false), check=false)
      set_attribute!(f, :inverse, f_inv)
      set_attribute!(f_inv, :inverse, f)
      f
    end
  )
  has_ancestor(W->any(WW->(WW === W), patches(C)), U) || error("patch not found")
  W = ambient_scheme(U)
  V = domain(iso)
  UV = codomain(iso)
  hV = complement_equation(UV)
  hU = complement_equation(U)
  WV = PrincipalOpenSubset(W, OO(W).([lifted_numerator(hU), lifted_numerator(hV)]))
  ident = SpecMor(UV, WV, hom(OO(WV), OO(UV), gens(OO(UV)), check=false), check=false)
  inv_ident = SpecMor(WV, UV, hom(OO(UV), OO(WV), gens(OO(WV)), check=false), check=false)
  new_iso =  compose(iso, ident)
  new_iso_inv = compose(inv_ident, inverse(iso))
  set_attribute!(new_iso, :inverse, new_iso_inv)
  set_attribute!(new_iso_inv, :inverse, new_iso)
  if any(WW->(WW===W), patches(C)) 
    return new_iso
  end
  return _flatten_open_subscheme(W, C, iso=new_iso)
end

function _flatten_open_subscheme(
    U::SimplifiedSpec, C::Covering;
    iso::AbsSpecMor=begin 
      UU = PrincipalOpenSubset(U, one(OO(U)))
      f = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false), check=false)
      f_inv = SpecMor(UU, U, hom(OO(U), OO(UU), gens(OO(UU)), check=false), check=false)
      set_attribute!(f, :inverse, f_inv)
      set_attribute!(f_inv, :inverse, f)
      f
    end
  )
  has_ancestor(W->any(WW->(WW === W), patches(C)), U) || error("patch not found")
  W = original(U)
  V = domain(iso)
  UV = codomain(iso)::PrincipalOpenSubset
  hV = complement_equation(UV)
  f, g = identification_maps(U)
  hVW = pullback(g)(hV)
  WV = PrincipalOpenSubset(W, hVW)
  ident = SpecMor(UV, WV, 
                  hom(OO(WV), OO(UV), 
                      OO(UV).(pullback(f).(gens(ambient_coordinate_ring(WV)))), 
                      check=false), 
                  check=false)
  inv_ident = SpecMor(WV, UV,
                      hom(OO(UV), OO(WV),
                          OO(WV).(pullback(g).(gens(ambient_coordinate_ring(UV)))),
                          check=false),
                      check=false)
  new_iso =  compose(iso, ident)
  new_iso_inv = compose(inv_ident, inverse(iso))
  set_attribute!(new_iso, :inverse, new_iso_inv)
  set_attribute!(new_iso_inv, :inverse, new_iso)
  if any(WW->(WW===W), patches(C)) 
    return new_iso
  end
  return _flatten_open_subscheme(W, C, iso=new_iso)
end

function _flatten_open_subscheme(
    U::AbsSpec, C::Covering;
    iso::AbsSpecMor=begin
      UU = PrincipalOpenSubset(U, one(OO(U)))
      f = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false), check=false)
      f_inv = SpecMor(UU, U, hom(OO(U), OO(UU), gens(OO(UU)), check=false), check=false)
      set_attribute!(f, :inverse, f_inv)
      set_attribute!(f_inv, :inverse, f)
      f
    end
  )
  return iso
end

function _flatten_open_subscheme(
    U::PrincipalOpenSubset, P::AbsSpec;
    iso::AbsSpecMor=begin
      UU = PrincipalOpenSubset(U, one(OO(U)))
      f = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false), check=false)
      f_inv = SpecMor(UU, U, hom(OO(U), OO(UU), gens(OO(UU)), check=false), check=false)
      set_attribute!(f, :inverse, f_inv)
      set_attribute!(f_inv, :inverse, f)
      f
    end
  )
  U === P && return iso

  has_ancestor(W->W===P, U) || error("ancestor not found")
  W = ambient_scheme(U)
  V = domain(iso)
  UV = codomain(iso)
  hV = complement_equation(UV)
  hU = complement_equation(U)
  WV = PrincipalOpenSubset(W, OO(W).([lifted_numerator(hU), lifted_numerator(hV)]))
  ident = SpecMor(UV, WV, hom(OO(WV), OO(UV), gens(OO(UV)), check=false), check=false)
  new_iso =  compose(iso, ident)
  new_iso_inv = compose(inverse(ident), inverse(iso))
  set_attribute!(new_iso, :inverse, new_iso_inv)
  set_attribute!(new_iso_inv, :inverse, new_iso)
  if W === P
    return new_iso
  end
  return _flatten_open_subscheme(W, P, iso=new_iso)
end

function _flatten_open_subscheme(
    U::SimplifiedSpec, P::AbsSpec;
    iso::AbsSpecMor=begin 
      UU = PrincipalOpenSubset(U, one(OO(U)))
      f = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false), check=false)
      f_inv = SpecMor(UU, U, hom(OO(U), OO(UU), gens(OO(UU)), check=false), check=false)
      set_attribute!(f, :inverse, f_inv)
      set_attribute!(f_inv, :inverse, f)
      f
    end
  )
  U === P && return iso

  has_ancestor(W->W===P, U) || error("ancestor not found")
  W = original(U)
  V = domain(iso)
  UV = codomain(iso)::PrincipalOpenSubset
  hV = complement_equation(UV)
  f, g = identification_maps(U)
  hVW = pullback(g)(hV)
  WV = PrincipalOpenSubset(W, hVW)
  ident = SpecMor(UV, WV, 
                  hom(OO(WV), OO(UV), 
                      OO(UV).(pullback(f).(gens(ambient_coordinate_ring(WV)))), 
                      check=false), 
                  check=false)
  new_iso =  compose(iso, ident)
  new_iso_inv = compose(inverse(ident), inverse(iso))
  set_attribute!(new_iso, :inverse, new_iso_inv)
  set_attribute!(new_iso_inv, :inverse, new_iso)
  if W === P
    return new_iso
  end
  return _flatten_open_subscheme(W, P, iso=new_iso)
end

function _flatten_open_subscheme(
    U::AbsSpec, P::AbsSpec;
    iso::AbsSpecMor=begin
      UU = PrincipalOpenSubset(U, one(OO(U)))
      f = SpecMor(U, UU, hom(OO(UU), OO(U), gens(OO(U)), check=false), check=false)
      f_inv = SpecMor(UU, U, hom(OO(U), OO(UU), gens(OO(UU)), check=false), check=false)
      set_attribute!(f, :inverse, f_inv)
      set_attribute!(f_inv, :inverse, f)
      f
    end
  )
  U === P || error("schemes have no valid relationship")
  return iso
end

########################################################################
# Lookup of common ancestors                                           #
#                                                                      #
# This returns the node in the natural tree structure which is closest #
# to the two input schemes.                                            #
########################################################################
function _have_common_ancestor(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
  U === V && return true, V
  if ambient_scheme(U) === V
    return true, V
  elseif ambient_scheme(V) === U
    return true, U
  elseif ambient_scheme(V) === ambient_scheme(U)
    return true, ambient_scheme(V)
  end
  return _have_common_ancestor(ambient_scheme(U), ambient_scheme(V))
end

function _have_common_ancestor(U::AbsSpec, V::PrincipalOpenSubset)
  return _have_common_ancestor(U, ambient_scheme(V))
end

function _have_common_ancestor(
    V::PrincipalOpenSubset,
    U::AbsSpec
  )
  return _have_common_ancestor(U, ambient_scheme(V))
end

function _have_common_ancestor(U::AbsSpec, V::AbsSpec)
  return U===V, U
end

function _have_common_ancestor(U::AbsSpec, V::SimplifiedSpec)
  return _have_common_ancestor(U, original(V))
end

function _have_common_ancestor(
    V::SimplifiedSpec,
    U::AbsSpec
  )
  return _have_common_ancestor(U, original(V))
end

function _have_common_ancestor(U::PrincipalOpenSubset, V::SimplifiedSpec)
  U === V && return true
  if ambient_scheme(U) === V
    return true, V
  elseif original(V) === U
    return true, U
  elseif original(V) === ambient_scheme(U)
    return true, original(V)
  end
  return _have_common_ancestor(ambient_scheme(U), original(V))
end

function _have_common_ancestor(
    V::SimplifiedSpec,
    U::PrincipalOpenSubset
  )
  return _have_common_ancestor(U, V)
end

function _have_common_ancestor(U::SimplifiedSpec, V::SimplifiedSpec)
  U === V && return true
  if original(U) === V
    return true, V
  elseif original(V) === U
    return true, U
  elseif original(V) === original(U)
    return true, original(V)
  end
  return _have_common_ancestor(original(U), original(V))
end
