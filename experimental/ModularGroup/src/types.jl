@attributes mutable struct ModularGroup
  s::PermGroupElem
  t::PermGroupElem
  r::PermGroupElem
  j::PermGroupElem
  
  function ModularGroup(s::PermGroupElem, t::PermGroupElem, r::PermGroupElem, j::PermGroupElem)
    if r != s^-1 * t^-1 * s || j != s^-1 * t^-1
      throw(ArgumentError("r and j are inconsistent with s and t"))
    end
    if !defines_coset_action_s_t(s, t)
      throw(ArgumentError("s and t do not describe the action of the generators S and T on the cosets of a finite-index subgroup of SL(2,Z)"))
    end
    return new(s, t, r, j)
  end
end
