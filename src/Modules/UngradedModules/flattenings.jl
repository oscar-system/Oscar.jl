########################################################################
# Rerouting of the module functionality via flattenings
########################################################################
function submodule_membership(::HasRingFlattening, a::FreeModElem, M::SubModuleOfFreeModule)
  F = ambient_free_module(M)
  Fb, iso_F = flatten(F)
  Mb, iso_M = flatten(M)
  return iso_F(a) in Mb
end

function _kernel(::HasRingFlattening, phi::FreeModuleHom{ModuleType, ModuleType}) where {ModuleType <: FreeMod}
  phi_b, iso_dom, iso_cod = flatten(phi)
  Kb, inc_Kb = kernel(phi_b)
  K, inc_K = sub(domain(phi), inverse(iso_dom).(gens(Kb)))
  #TODO: Cache the flattenings of K and inc_K.
  return K, inc_K
end
