function generalized_cyclic_torus_cover(n::Int64, d::Int64, vslits::Vector, hslits::Vector)
  return GAP.Globals.GeneralizedCyclicTorusCover(n, d, GapObj(vslits), GapObj(hslits))
end

function comb_origami(n::Int64, x::Int64, y::Int64)
  return GAP.Globals.CombOrigami(n, x, y)
end

function cyclic_torus_cover_origamiS(n::Int64, d::Int64, v::Vector)
  return GAP.Globals.CyclicTorusCoverOrigamiS(n, d, GapObj(v))
end

function cyclic_torus_cover_origamiL(n::Int64, d::Int64, v::Vector)
  return GAP.Globals.CyclicTorusCoverOrigamiL(n, d, GapObj(v))
end

function base_change_l_to_s(n::Int64)
  return GAP.Globals.BaseChangeLToS(n)
end

function translation_group_on_homology_of_tn(n::Int64)
  return GAP.Globals.TranslationGroupOnHomologyOfTn(n)
end

function action_of_t_on_homology_of_tn(n::Int64)
  return GAP.Globals.ActionOfTOnHomologyOfTn(n)
end

function action_of_s_on_homology_of_tn(n::Int64)
  return GAP.Globals.ActionOfSOnHomologyOfTn(n)
end

function action_of_matrix_on_homology_of_tn(n::Int64, A::Matrix)
  return GAP.Globals.ActionOfMatrixOnHomologyOfTn(n, GapObj(A))
end
