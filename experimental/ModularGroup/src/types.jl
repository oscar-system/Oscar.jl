@attributes mutable struct ModularGroup <: Group
  s::PermGroupElem
  t::PermGroupElem
  
  function ModularGroup(s::PermGroupElem, t::PermGroupElem, check::Bool)
    if check
      @req defines_coset_action_s_t(s, t) "s and t do not describe the action of the generators S and T on the cosets of a finite-index subgroup of SL(2,Z)"
    end
    return new(s, t)
  end
end

struct ModularGroupElem <: GroupElem
  parent::ModularGroup
  mat::ZZMatrix
end
