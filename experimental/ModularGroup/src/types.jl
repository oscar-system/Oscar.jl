@attributes mutable struct ModularGroup
  s::PermGroupElem
  t::PermGroupElem
  r::PermGroupElem
  j::PermGroupElem
  ModularGroup(s::PermGroupElem, t::PermGroupElem, r::PermGroupElem, j::PermGroupElem) = new(s, t, r, j)
end
