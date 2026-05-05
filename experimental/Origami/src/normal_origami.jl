function normal_stored_origami(h::PermGroupElem, v::PermGroupElem, g::Oscar.GAPGroup)
  return Normal_stored_origami(
    GAP.Globals.NormalStoredOrigami(GapObj(h), GapObj(v), g), h, v, g
  )
end
function as_permutation_pepresentation(o::Normal_stored_origami)
  return GAP.Globals.AsPermutationRepresentation(GapObj(o))
end

function all_normal_origamis_by_degree(n::Int)
  return GAP.Globals.AllNormalOrigamisByDegree(n)
end

function all_normal_origamis_from_group(g::GAPGroup)
  return GAP.Globals.AllNormalOrigamisFromGroup(GapObj(g))
end
