# We implement the SubModuleOfFreeModule-level of the module functionality 
# interface for modules over MPolyQuos.

########################################################################
# The essential three functions:                                       #
########################################################################
@attr function kernel(
    f::FreeModuleHom{DomainType, CodomainType}
  ) where {
           DomainType<:FreeMod{<:MPolyQuoElem},
           CodomainType<:FreeMod{<:MPolyQuoElem}
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
  KK, inc2 = sub(domain(f), tr.(gens(K)))
  return KK, inc2
end

function Base.in(a::FreeModElem{T}, M::SubModuleOfFreeModule{T}) where {T<:MPolyQuoElem}
  return _lifting_map(parent(a))(a) in _poly_module(M)
end

function coordinates(
    v::FreeModElem{T}, 
    M::SubModuleOfFreeModule{T}
  ) where {T<:MPolyQuoElem}
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
@attr function kernel(
    f::FreeModuleHom{DomainType, CodomainType}
  ) where {
           DomainType<:FreeMod{<:MPolyQuoElem},
           CodomainType<:SubQuo{<:MPolyQuoElem}
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
  KK, inc2 = sub(domain(f), tr.(gens(K)))
  return KK, inc2
end

function coordinates(
    v::FreeModElem{T}, 
    M::SubModuleOfFreeModule{T}, 
    task::Symbol
  ) where {T<:MPolyQuoElem}
  return coordinates(v, M)
end

#function free_resolution(M::SubQuo{T}) where {T<:MPolyQuoElem}
#  R = base_ring(M)
#  p = presentation(M)
#  K, inc = kernel(map(p, 1))
#  i = 1
#  while !iszero(K)
#    F = FreeMod(R, ngens(K))
#    phi = hom(F, p[i], inc.(gens(K)))
#    p = Hecke.ChainComplex(ModuleFP, pushfirst!(ModuleFPHom[map(p, i) for i in collect(range(p))[1:end-1]], phi), check=false, seed = -2)
#    i = i+1
#    K, inc = kernel(phi)
#  end
#  #end_map = hom(FreeMod(R, 0), K, elem_type(K)[])
#  p = Hecke.ChainComplex(ModuleFP, vcat(ModuleFPHom[inc], ModuleFPHom[map(p, i) for i in collect(range(p))[1:end-1]]), check=false, seed = -2)
#  return p
#end


########################################################################
# Auxiliary helping functions to allow for the above                   #
########################################################################
#
### For a free module F = R^r over R = P/I, this returns a lifting map 
# to the module P^r/I*P^r. Note that this is an unnatural map since 
# the latter is an R-module only by accident.
@attr function _lifting_iso(F::FreeMod{T}) where {T<:MPolyQuoElem}
  M = _as_poly_module(F)
  function my_lift(v::FreeModElem{T}) where {T<:MPolyQuoElem}
    parent(v) === F || error("element does not have the right parent")
    w = elem_type(M)[lift(a)*M[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(M)
    return sum(w)
  end
  return my_lift
end

### For a free module F = R^r over R = P/I, this returns a lifting map 
# to the module P^r. Note that this is not a homomorphism of modules. 
@attr function _lifting_map(F::FreeMod{T}) where {T<:MPolyQuoElem}
  FP = _poly_module(F)
  function my_lift(v::FreeModElem{T}) where {T<:MPolyQuoElem}
    parent(v) === F || error("element does not have the right parent")
    w = [lift(a)*FP[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(FP)
    return sum(w)
  end
  return my_lift
end

### To a free module over R = P/I, return the free module over R 
# in the same number of generators
@attr function _poly_module(F::FreeMod{T}) where {T<:MPolyQuoElem}
  R = base_ring(F)
  P = base_ring(R) # the polynomial ring
  r = rank(F)
  FP = FreeMod(P, r) 
  return FP
end

### Return the canonical projection FP -> F from the P-module FP to the 
# R-module F.
@attr function _poly_module_restriction(F::FreeMod{T}) where {T<:MPolyQuoElem}
  R = base_ring(F)
  P = base_ring(R)
  FP = _poly_module(F)
  return hom(FP, F, gens(F), x->R(x))
end

### Return the same module, but as a SubQuo over the polynomial ring
@attr function _as_poly_module(F::FreeMod{T}) where {T<:MPolyQuoElem}
  R = base_ring(F)
  P = base_ring(R)
  I = modulus(R)
  FP = FreeMod(P, rank(F))
  IFP, inc = I*FP
  M, p = quo(FP, IFP)
  return M
end

### Return an isomorphism with _as_poly_module(F)
@attr function _iso_with_poly_module(F::FreeMod{T}) where {T<:MPolyQuoElem}
  M = _as_poly_module(F)
  return hom(M, F, gens(F), x->(base_ring(F)(x)))
end

### Return the preimage of M under the canonical projection P^r -> R^r 
# for R^r the ambient_free_module of M.
@attr function _poly_module(M::SubModuleOfFreeModule{T}) where {T<:MPolyQuoElem}
  F = ambient_free_module(M) 
  FP = _poly_module(F)
  v = [_lifting_map(F)(g) for g in gens(M)] 
  w = [f*e for e in gens(FP) for f in gens(modulus(base_ring(M)))]
  MP = SubModuleOfFreeModule(FP, vcat(v, w))
  return MP
end

@attr function _as_poly_module(M::SubQuo{T}) where {T<:MPolyQuoElem}
  F = ambient_free_module(M) 
  FP = _poly_module(F)
  v = [_lifting_map(F)(g) for g in ambient_representatives_generators(M)] 
  w = [f*e for e in gens(FP) for f in gens(modulus(base_ring(M)))]
  w_ext = vcat(w, [_lifting_map(F)(g) for g in relations(M)])
  MP = SubQuo(FP, v, w_ext)
  return MP
end

@attr function _iso_with_poly_module(F::SubQuo{T}) where {T<:MPolyQuoElem}
  M = _as_poly_module(F)
  return hom(M, F, gens(F), x->(base_ring(F)(x)))
end

@attr function _lifting_iso(F::SubQuo{T}) where {T<:MPolyQuoElem}
  M = _as_poly_module(F)
  function my_lift(v::SubQuoElem{T}) where {T<:MPolyQuoElem}
    parent(v) === F || error("element does not have the right parent")
    w = elem_type(M)[lift(a)*M[i] for (i, a) in coordinates(v)]
    iszero(length(w)) && return zero(M)
    return sum(w)
  end
  return my_lift
end
