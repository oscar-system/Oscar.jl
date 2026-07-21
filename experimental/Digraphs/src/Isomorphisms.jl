function is_isomorphic(D1::Digraph, D2::Digraph)
  return DigraphWrap.IsIsomorphicDigraph(GapObj(D1), GapObj(D2))::Bool
end

function automorphism_group(D::Digraph)
  return PermGroup(DigraphWrap.AutomorphismGroup(GapObj(D)))
end

function bliss_automorphism_group(D::Digraph)
  return PermGroup(DigraphWrap.BlissAutomorphismGroup(GapObj(D)))
end

function bliss_canonical_labelling(D::Digraph)
  return DigraphWrap.BlissCanonicalLabelling(GapObj(D))
end


