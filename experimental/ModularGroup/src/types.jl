mutable struct ModularGroup
  s::PermGroupElem
  t::PermGroupElem
  r::PermGroupElem
  j::PermGroupElem
  perm_group_cache::Union{Nothing, PermGroup}
  coset_action_hom_cache::Union{Nothing, GAPGroupHomomorphism{FPGroup, PermGroup}}
end
