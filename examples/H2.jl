module H2_G_QmodZ_mod
using Oscar

function schur_cover(G::Oscar.GAPGroup)
  f = GAP.Globals.EpimorphismSchurCover(G.X)
  k = GAP.Globals.Source(f)
  S = FPGroup(k)
  return S, GAPGroupHomomorphism(S, G, f)
end
  

#from https://sites.google.com/view/andre-macedo/code?pli=1
"""
Should compute
    H^2(G, Q/Z)
for "any" Gap-group G. For abelian (GrpAbFinGen), there is a more
direct implementation in Hecke.
"""
function H2_G_QmodZ(G::Oscar.GAPGroup)
  _, f = schur_cover(G)
  return kernel(f)[1]
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
  _, f = schur_cover(G)
  k = kernel(f)[1]
  im_gen = elem_type(k)[]
  for D = U
    gp = intersect(k, derived_subgroup(preimage(f, D)[1])[1])[1]
    for g = gens(gp)
      push!(im_gen, g)
    end
  end
  #this works:
#  return FPGroup(k.X/GAP.Globals.Subgroup(k.X, GAP.GapObj(im_gen)))
  #this doesn't:
  #= example input
  A = abelian_group(PermGroup, [2,2])
  B = [sub(A, [A[1]])[1], sub(A, [A[2]])[1]]
  H2_G_QmodZ_kern_restriction(A, B)
  =#

  return quo(k, im_gen)
end  


end #module H2_G_QmodZ_mod

