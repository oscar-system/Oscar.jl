# We implement the SubModuleOfFreeModule-level of the module functionality 
# interface for modules over MPolyQuos.

########################################################################
# The essential three functions:                                       #
########################################################################
@attr Any function kernel(
    f::FreeModuleHom{DomainType, CodomainType}
  ) where {
           DomainType<:FreeMod{<:MPolyQuoRingElem},
           CodomainType<:FreeMod{<:MPolyQuoRingElem}
          }
  R = base_ring(codomain(f))
  P = base_ring(R)
  F = _poly_module(domain(f))
  M = _as_poly_module(codomain(f))
  id = _iso_with_poly_module(codomain(f))
  # Why does img_gens(f) return a list of SubQuoElems???
  phi = _lifting_iso(codomain(f))
  g = hom(F, M, phi.(f.(gens(domain(f)))))
  K, inc = kernel(g)
  tr =  compose(inc, _poly_module_restriction(domain(f)))
  KK, inc2 = sub(domain(f), unique!(filter!(!iszero, tr.(gens(K)))))
  return KK, inc2
end

function Base.in(a::FreeModElem{T}, M::SubModuleOfFreeModule{T}) where {T<:MPolyQuoRingElem}
  return _lifting_map(parent(a))(a) in _poly_module(M)
end

function coordinates(
    v::FreeModElem{T}, 
    M::SubModuleOfFreeModule{T}
  ) where {T<:MPolyQuoRingElem}
  w = _lifting_map(parent(v))(v)
  MP = _poly_module(M)
  c = coordinates(w, MP)
  R = base_ring(M)
  entries = [(i, R(a)) for (i, a) in c if i <= ngens(M)]
  return sparse_row(R, entries)
end

########################################################################
# Methods which should not be necessary, but the stuff doesn't work,   #
# unless we implement them.                                            #
########################################################################

function coordinates(
    v::FreeModElem{T}, 
    M::SubModuleOfFreeModule{T}, 
    task::Symbol
  ) where {T<:MPolyQuoRingElem}
  return coordinates(v, M)
end

########################################################################
# Auxiliary helping functions to allow for the above                   #
########################################################################
#
### For a free module F = R^r over R = P/I, this returns a lifting map 
# to the module P^r/I*P^r. Note that this is an unnatural map since 
# the latter is an R-module only by accident.
@attr Any function _lifting_iso(F::FreeMod{T}) where {T<:MPolyQuoRingElem}
  M = _as_poly_module(F)
  function my_lift(v::FreeModElem{T}) where {T<:MPolyQuoRingElem}
    parent(v) === F || error("element does not have the right parent")
    w = elem_type(M)[lift(a)*M[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(M)
    return sum(w)
  end
  return my_lift
end

### For a free module F = R^r over R = P/I, this returns a lifting map 
# to the module P^r. Note that this is not a homomorphism of modules. 
@attr Any function _lifting_map(F::FreeMod{T}) where {T<:MPolyQuoRingElem}
  FP = _poly_module(F)
  function my_lift(v::FreeModElem{T}) where {T<:MPolyQuoRingElem}
    parent(v) === F || error("element does not have the right parent")
    w = [lift(a)*FP[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(FP)
    return sum(w)
  end
  return my_lift
end

### To a free module over R = P/I, return the free module over R 
# in the same number of generators
@attr Any function _poly_module(F::FreeMod{T}) where {T<:MPolyQuoRingElem}
  R = base_ring(F)
  P = base_ring(R) # the polynomial ring
  r = rank(F)
  FP = FreeMod(P, r) 
  return FP
end

### Return the canonical projection FP -> F from the P-module FP to the 
# R-module F.
@attr Any function _poly_module_restriction(F::FreeMod{T}) where {T<:MPolyQuoRingElem}
  R = base_ring(F)
  P = base_ring(R)
  FP = _poly_module(F)
  return hom(FP, F, gens(F), R)
end

### Return the same module, but as a SubquoModule over the polynomial ring
@attr Any function _as_poly_module(F::FreeMod{T}) where {T<:MPolyQuoRingElem}
  R = base_ring(F)
  P = base_ring(R)
  I = modulus(R)
  FP = FreeMod(P, rank(F))
  IFP, inc = I*FP
  M, p = quo(FP, IFP)
  # Manually set a groebner basis
  # This is brittle! Please improve!
  J = ideal(P, gens(groebner_basis(I))) # a groebner basis has already been computed.
  JFP, _ = J*FP
  G, _ = quo(FP, JFP)
  gb = G.quo.gens
  ord = default_ordering(M.quo)
  gb.ordering = ord
  singular_assure(gb)
  gb.isGB = true
  gb.S.isGB = true
  M.quo.groebner_basis[ord] = gb
  return M
end

### Return an isomorphism with _as_poly_module(F)
@attr Any function _iso_with_poly_module(F::FreeMod{T}) where {T<:MPolyQuoRingElem}
  M = _as_poly_module(F)
  return hom(M, F, gens(F), base_ring(F))
end

### Return the preimage of M under the canonical projection P^r -> R^r 
# for R^r the ambient_free_module of M.
@attr Any function _poly_module(M::SubModuleOfFreeModule{T}) where {T<:MPolyQuoRingElem}
  F = ambient_free_module(M) 
  FP = _poly_module(F)
  v = elem_type(FP)[_lifting_map(F)(g) for g in gens(M)] 
  w = elem_type(FP)[f*e for e in gens(FP) for f in gens(modulus(base_ring(M)))]
  MP = SubModuleOfFreeModule(FP, vcat(v, w))
  return MP
end

@attr Any function _as_poly_module(M::SubquoModule{T}) where {T<:MPolyQuoRingElem}
  F = ambient_free_module(M) 
  FP = _poly_module(F)
  v = [_lifting_map(F)(g) for g in ambient_representatives_generators(M)] 
  w = [f*e for e in gens(FP) for f in gens(modulus(base_ring(M)))]
  w_ext = vcat(w, elem_type(FP)[_lifting_map(F)(g) for g in relations(M)])
  MP = SubquoModule(FP, v, w_ext)
  return MP
end

@attr Any function _iso_with_poly_module(F::SubquoModule{T}) where {T<:MPolyQuoRingElem}
  M = _as_poly_module(F)
  return hom(M, F, gens(F), base_ring(F))
end

@attr Any function _lifting_iso(F::SubquoModule{T}) where {T<:MPolyQuoRingElem}
  M = _as_poly_module(F)
  function my_lift(v::SubquoModuleElem{T}) where {T<:MPolyQuoRingElem}
    parent(v) === F || error("element does not have the right parent")
    w = elem_type(M)[lift(a)*M[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(M)
    return sum(w)
  end
  return my_lift
end
