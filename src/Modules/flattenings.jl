### Flattenings of (graded) modules

function flatten(F::ModuleFP)
  return flatten(base_ring(F))(F)
end

function (flat_map::RingFlattening)(F::ModuleFP{T}) where {T <: MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                                                  MPolyQuoLocRingElem, MPolyLocRingElem}}}
  if !haskey(flat_counterparts(flat_map), F)
    F_flat, iso = change_base_ring(flat_map, F)
    iso_inv = hom(F_flat, F, gens(F), inverse(flat_map))
    set_attribute!(iso, :inverse, iso_inv)
    set_attribute!(iso_inv, :inverse, iso)
    flat_counterparts(flat_map)[F] = (F_flat, iso, iso_inv)
  end
    
  F_flat, iso, iso_inv = flat_counterparts(flat_map)[F]::Tuple{<:ModuleFP, <:ModuleFPHom, <:ModuleFPHom}
  return F_flat, iso, iso_inv
end

function flatten(
    a::FreeModElem{T}
  ) where {T <: MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                      MPolyQuoLocRingElem, MPolyLocRingElem}}}
  F = parent(a)
  R = base_ring(F)
  flat_map = flatten(R)
  return flat_map(a)
end

function (flat::RingFlattening)(a::FreeModElem)
  if !haskey(flat_counterparts(flat), a)
    F = parent(a)
    F_flat, iso, iso_inv = flat(F)
    a_flat = iso(a)
    flat_counterparts(flat)[a] = a_flat
  end
  return flat_counterparts(flat)[a]::FreeModElem
end


function Base.in(
    a::FreeModElem{T}, M::SubquoModule{T}
  ) where {T <: MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                      MPolyQuoLocRingElem, MPolyLocRingElem}}}
  flat = flatten(base_ring(parent(a)))
  return flat(a) in flat(M)[1]
end

function coordinates(
    a::FreeModElem{T}, M::SubquoModule{T}
  ) where {T <: MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                      MPolyQuoLocRingElem, MPolyLocRingElem}}}

  flat = flatten(base_ring(parent(a)))
  c = coordinates(flat(a), flat(M)[1])
  return map_entries(inverse(flat), c)
end

function free_resolution(
    M::SubquoModule{T}
  ) where {T <: MPolyRingElem{<:Union{MPolyRingElem, MPolyQuoRingElem, 
                                      MPolyQuoLocRingElem, MPolyLocRingElem}}}
  flat = flatten(base_ring(M))
  M_flat, iso_M, iso_M_inv = flat(M)
  comp = free_resolution(M_flat) # assuming that this is potentially cached
  if !haskey(flat_counterparts(flat), comp)
    res_obj = ModuleFP[]
    isos = ModuleFPHom[]
    push!(res_obj, M)
    push!(isos, iso_M_inv)
    res_maps = Map[]
    for i in 0:first(chain_range(comp))
      F_flat = comp[i]
      @assert F_flat === domain(map(comp, i))
      @assert domain(last(isos)) === codomain(map(comp, i))
      F, iso_F_inv = _change_base_ring_and_preserve_gradings(inverse(flat), F_flat)
      iso_F = hom(F, F_flat, gens(F_flat), flat)
      push!(res_maps, hom(F, last(res_obj), last(isos).(map(comp, i).(iso_F.(gens(F))))))
      push!(res_obj, F)
      push!(isos, iso_F_inv)
    end
    comp_up = ComplexOfMorphisms(ModuleFP, reverse(res_maps), typ=:chain, seed=-1, check=false)
    comp_up.complete = true
    result = FreeResolution(comp_up)
    flat_counterparts(flat)[comp] = result
  end
  return flat_counterparts(flat)[comp]::FreeResolution
end

# TODO: We need this special routine, because we can not write a generic 
# base change method for morphisms of graded rings. See pr #2677
function _change_base_ring_and_preserve_gradings(phi::Any, F::FreeMod)
  R = base_ring(F)
  S = parent(phi(zero(R)))
  FS = (is_graded(F) ? graded_free_module(S, degree.(gens(F))) : FreeMod(S, ngens(F)))
  FS.S = F.S
  return FS, hom(F, FS, gens(FS), phi)
end

function _change_base_ring_and_preserve_gradings(phi::Any, M::SubquoModule)
  R = base_ring(M)
  S = parent(phi(zero(R)))
  F = ambient_free_module(M)
  FF, iso_F = _change_base_ring_and_preserve_gradings(phi, F)
  MM = SubquoModule(FF, iso_F.(ambient_representatives_generators(M)), iso_F.(relations(M)))
  return MM, hom(M, MM, gens(MM), phi)
end

function _change_base_ring_and_preserve_gradings(
    phi::Any, f::ModuleFPHom;
    domain_change::Map = _change_base_ring_and_preserve_gradings(phi, domain(f))[2],
    codomain_change::Map = _change_base_ring_and_preserve_gradings(phi, codomain(f))[2]
  )
  D = codomain(domain_change)
  C = codomain(codomain_change)
  result = hom(D, C, codomain_change.(f.(gens(domain(f)))))
  return result
end

