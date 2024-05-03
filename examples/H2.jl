module H2_G_QmodZ_mod
using Oscar

function schur_cover(G::Oscar.GAPGroup)
  f = GAP.Globals.EpimorphismSchurCover(GapObj(G))
  k = GAP.Globals.Source(f)
  S = FPGroup(k)
  return S, GAPGroupHomomorphism(S, G, f)
end
  

#from https://sites.google.com/view/andre-macedo/code?pli=1
"""
Should compute
    H^2(G, Q/Z)
for "any" Gap-group G. For abelian (FinGenAbGroup), there is a more
direct implementation in Hecke.
"""
function H2_G_QmodZ(G::Oscar.GAPGroup)
  _, f = schur_cover(G)
  return kernel(f)[1]
end

"""
Should compute
    Kern(H^2(G, Q/Z) -> sum_(u in U) H^2(u, Q/Z))
Where the "->" is the restriction from H^2(G) -> H^2(U).

Useful(?) is U is a family of decomposition groups.
"""
function H2_G_QmodZ_kern_restriction(G::T, U::Vector{T}) where T <: Oscar.GAPGroup
  _, f = schur_cover(G)
  k = kernel(f)[1]
  im_gen = elem_type(k)[]
  for D = U
    gp = intersect(k, derived_subgroup(preimage(f, D)[1])[1])[1]
    for g = gens(gp)
      push!(im_gen, g)
    end
  end
  return quo(k, im_gen)
end  


end #module H2_G_QmodZ_mod

