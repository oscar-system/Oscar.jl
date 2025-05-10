### Flattenings of (graded) modules over polynomial rings and quotients of polynomial rings

function flatten(F::SparseFPModule)
  return flatten(base_ring(F))(F)
end

const FlattableRingElemType = Union{<:MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                                            MPolyQuoLocRingElem, MPolyLocRingElem}},
                                    <:MPolyQuoRingElem{<:MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                                                               MPolyQuoLocRingElem, MPolyLocRingElem}}}
                                   }

function (flat_map::RingFlattening)(F::FreeMod{T}) where {T <: FlattableRingElemType}
  return get!(flat_counterparts(flat_map), F) do
    F_flat, iso = _change_base_ring_and_preserve_gradings(flat_map, F)
    iso_inv = hom(F_flat, F, gens(F), inverse(flat_map); check=false)
    set_attribute!(iso, :inverse, iso_inv)
    set_attribute!(iso_inv, :inverse, iso)
    return F_flat, iso, iso_inv
  end::Tuple{<:SparseFPModule, <:SparseFPModuleHom, <:SparseFPModuleHom}
end

function (flat_map::RingFlattening)(M::SubquoModule{T}) where {T <: FlattableRingElemType}
  return get!(flat_counterparts(flat_map), M) do
    F = ambient_free_module(M)
    F_flat, iso_F, iso_inv_F = flat_map(F)
    M_flat = SubquoModule(F_flat, 
                          iso_F.(ambient_representatives_generators(M)), 
                          iso_F.(relations(M))
                         )
    iso = hom(M, M_flat, gens(M_flat), flat_map; check=false)
    iso_inv = hom(M_flat, M, gens(M), inverse(flat_map); check=false)
    set_attribute!(iso, :inverse, iso_inv)
    set_attribute!(iso_inv, :inverse, iso)
    return M_flat, iso, iso_inv
  end::Tuple{<:SparseFPModule, <:SparseFPModuleHom, <:SparseFPModuleHom}
end

function flatten(
    a::FreeModElem{T}
  ) where {T <: FlattableRingElemType}
  F = parent(a)
  R = base_ring(F)
  flat_map = flatten(R)
  return flat_map(a)
end

function (flat::RingFlattening)(a::FreeModElem)
  return get!(flat_counterparts(flat), a) do
    F = parent(a)
    F_flat, iso, iso_inv = flat(F)
    return iso(a)
  end::FreeModElem
end

### TODO: The following two functions should not be necessary if the module code 
# was successfully referring everything to the SubModuleOfFreeModule-layer. 
# Once these layers are finally documented, this should be cleaned up
function Base.in(
    a::FreeModElem{T}, M::SubquoModule{T}
  ) where {T <: FlattableRingElemType}
  flat = flatten(base_ring(parent(a)))
  return flat(a) in flat(M)[1]
end

function coordinates(
    a::FreeModElem{T}, M::SubquoModule{T}
  ) where {T <: FlattableRingElemType}

  flat = flatten(base_ring(parent(a)))
  c = coordinates(flat(a), flat(M)[1])
  is_zero(c) && return sparse_row(base_ring(a))
  return map_entries(inverse(flat), c)
end

function free_resolution(
    M::SubquoModule{T}
  ) where {T <: FlattableRingElemType}
  flat = flatten(base_ring(M))
  M_flat, iso_M, iso_M_inv = flat(M)
  comp = free_resolution(M_flat) # assuming that this is potentially cached
  return get!(flat_counterparts(flat), comp) do
    res_obj = SparseFPModule[]
    isos = SparseFPModuleHom[]
    push!(res_obj, M)
    push!(isos, iso_M_inv)
    res_maps = Map[]
    for i in 0:first(chain_range(comp))
      F_flat = comp[i]
      @assert F_flat === domain(map(comp, i))
      @assert domain(last(isos)) === codomain(map(comp, i))
      F, iso_F_inv = _change_base_ring_and_preserve_gradings(inverse(flat), F_flat)
      iso_F = hom(F, F_flat, gens(F_flat), flat; check=false)
      push!(res_maps, hom(F, last(res_obj), last(isos).(map(comp, i).(iso_F.(gens(F)))); check=false))
      push!(res_obj, F)
      push!(isos, iso_F_inv)
    end
    comp_up = ComplexOfMorphisms(SparseFPModule, reverse(res_maps), typ=:chain, seed=-1, check=false)
    comp_up.complete = true
    return FreeResolution(comp_up)
  end::FreeResolution
end

function (phi::RingFlattening)(f::SparseFPModuleHom)
  return get!(flat_counterparts(phi), f) do
    dom = domain(f)
    dom_b, iso_dom, iso_inv_dom = phi(dom)
    cod = codomain(f)
    cod_b, iso_cod, iso_inv_cod = phi(cod)
    fb = hom(dom_b, cod_b, iso_cod.(f.(gens(dom))); check=false)
    return fb
  end::SparseFPModuleHom
end

# TODO: We need this special routine, because we can not write a generic 
# base change method for morphisms of graded rings. See pr #2677
function _change_base_ring_and_preserve_gradings(phi::Any, F::FreeMod)
  R = base_ring(F)
  S = parent(phi(zero(R)))
  FS = (is_graded(F) ? graded_free_module(S, [degree(g; check=false) for g in gens(F)]) : FreeMod(S, ngens(F)))
  FS.S = F.S
  return FS, hom(F, FS, gens(FS), phi; check=false)
end

function _change_base_ring_and_preserve_gradings(
    phi::Any, M::SubquoModule;
    ambient_base_change::FreeModuleHom=begin
      F = ambient_free_module(M)
      FF, iso_F = _change_base_ring_and_preserve_gradings(phi, F)
      iso_F
    end
  )
  FF = codomain(ambient_base_change)
  MM = SubquoModule(FF, ambient_base_change.(ambient_representatives_generators(M)), 
                    ambient_base_change.(relations(M)))
  return MM, hom(M, MM, gens(MM), phi; check=false)
end

function _change_base_ring_and_preserve_gradings(
    phi::Any, f::SparseFPModuleHom;
    domain_change::Map = _change_base_ring_and_preserve_gradings(phi, domain(f))[2],
    codomain_change::Map = _change_base_ring_and_preserve_gradings(phi, codomain(f))[2]
  )
  D = codomain(domain_change)
  C = codomain(codomain_change)
  result = hom(D, C, codomain_change.(f.(gens(domain(f)))); check=false)
  return result
end

# SubModuleOfFreeModule needs its own treatment as they differ from the user facing 
# modules in that they don't have their own elements and they dont have their own 
# homomorphisms.
function (flat_map::RingFlattening)(
                                    I::SubModuleOfFreeModule{T}
                                   ) where {T <: FlattableRingElemType}
  return get!(flat_counterparts(flat_map), I) do
    F = ambient_free_module(I)
    R = base_ring(I)
    flat_map = flatten(R)
    Fb, iso_F = flat_map(F)
    return _change_base_ring_and_preserve_gradings(flat_map, I; ambient_base_change=iso_F)
  end::SubModuleOfFreeModule
end

function _change_base_ring_and_preserve_gradings(
    phi::Any, M::SubModuleOfFreeModule;
    ambient_base_change::FreeModuleHom=begin
      F = ambient_free_module(M)
      FF, iso_F = _change_base_ring_and_preserve_gradings(phi, F)
      iso_F
    end
  )
  FF = codomain(ambient_base_change)
  MM = SubModuleOfFreeModule(FF, ambient_base_change.(gens(M)))
  # These modules dont have their own homomorphisms. 
  # Hence the return signature differs.
  return MM
end

### The magic three functions for the SubModuleOfFreeModule layer for flattenings.
function Base.in(
    a::FreeModElem{T}, M::SubModuleOfFreeModule{T}
  ) where {T <: FlattableRingElemType}
  flat = flatten(base_ring(parent(a)))
  return flat(a) in flat(M)
end

function coordinates(
    a::FreeModElem{T}, M::SubModuleOfFreeModule{T}
  ) where {T <: FlattableRingElemType}
  flat = flatten(base_ring(parent(a)))
  c = coordinates(flat(a), flat(M))
  is_zero(c) && return sparse_row(base_ring(a))
  return map_entries(inverse(flat), c)
end

function kernel(
    phi::FreeModuleHom{
                       <:FreeMod{<:FlattableRingElemType},
                       <:FreeMod{<:FlattableRingElemType},
                     Nothing
                    }
  )
  R = base_ring(domain(phi))
  flat_map = flatten(R)
  dom_flat, iso_dom, iso_dom_inv = flat_map(domain(phi))
  cod_flat, iso_cod, iso_cod_inv = flat_map(codomain(phi))
  phi_b = flat_map(phi)
  Kb, inc_Kb = kernel(phi_b)
  K, inc_K = sub(domain(phi), iso_dom_inv.(ambient_representatives_generators(Kb)))
  iso_K = hom(K, Kb, gens(Kb), flat_map; check=false)
  iso_inv_K = hom(Kb, K, gens(K), inverse(flat_map); check=false)
  flat_counterparts(flat_map)[K] = Kb, iso_K, iso_inv_K
  flat_counterparts(flat_map)[inc_K] = inc_Kb
  return K, inc_K
end

