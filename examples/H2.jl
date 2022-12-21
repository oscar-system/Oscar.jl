module H2_G_QmodZ_mod
using Oscar

#from https://sites.google.com/view/andre-macedo/code?pli=1
"""
Should compute
    H^2(G, Q/Z)
for "any" Gap-group G. For abelian (GrpAbFinGen), there is a more
direct implementation in Hecke.
"""
function H2_G_QmodZ(G::Oscar.GAPGroup)
  f = GAP.Globals.EpimorphismSchurCover(G.X)
  k = GAP.Globals.Kernel(f)
  return FPGroup(k)
end

"""
Should compute
    Kern(H^2(G, Q/Z) -> sum_(u in U) H^2(u, Q/Z))
Where the "->" is the restriction from H^2(G) -> H^2(U).

Useful(?) is U is a family of decomposition groups.
"""
function H2_G_QmodZ_kern_restriction(G::T, U::Vector{T}) where T <: Oscar.GAPGroup
  f = GAP.Globals.EpimorphismSchurCover(G.X)
  k = GAP.Globals.Kernel(f)
  im_gen = []
  for D = U
    gp = GAP.Globals.Intersection(k,GAP.Globals.DerivedSubgroup(GAP.Globals.PreImagesSet(f, D.X)))
    for g = GAP.Globals.GeneratorsOfGroup(gp)
      push!(im_gen, g)
    end
  end
  return FPGroup(k/GAP.Globals.Subgroup(k, GAP.GapObj(im_gen)))
end  


end #module H2_G_QmodZ_mod

