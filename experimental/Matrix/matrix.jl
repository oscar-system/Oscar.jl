module MatrixGroups

using GAP
using Oscar

export wrap_for_gap
  
function __init__()
    GAP.Globals.Reread(GAP.GapObj(joinpath(Oscar.oscardir, "experimental", "Matrix", "matrix.g")))
end

function _wrap_for_gap(m::MatrixElem)
    return GAP.Globals.MakeJuliaMatrixRep(m)
end

function MatrixGroup(matrices::Vector{<:MatrixElem{T}}) where T <: Union{fmpz, fmpq, nf_elem}
       @assert !isempty(matrices)
    
       # One should probably check whether all matrices are n by n (and invertible
       # and such ...)

       K = base_ring(matrices[1])
       if K isa FlintIntegerRing
          K = QQ
       end
       n = nrows(matrices[1])

       Fq, matrices_Fq, OtoFq = Oscar.good_reduction(matrices, 2)

       G = Oscar.MatrixGroup(n, Fq, matrices_Fq)
       N = order(G)
       if !isdivisible_by(Hecke._minkowski_multiple(K, n), N)
          error("Group is not finite")
       end

       #G_to_fin_pres = GAP.Globals.IsomorphismFpGroupByGenerators(G.X, GapObj([ g.X for g in gens(G) ]))
       #F = GAP.Globals.Range(G_to_fin_pres)
       #rels = GAP.Globals.RelatorsOfFpGroup(F)

       #gens_and_invsF = [ g for g in GAP.Globals.FreeGeneratorsOfFpGroup(F) ]
       #append!(gens_and_invsF, [ inv(g) for g in GAP.Globals.FreeGeneratorsOfFpGroup(F) ])
       #matrices_and_invs = copy(matrices)
       #append!(matrices_and_invs, [ inv(M) for M in matrices ])
       #for i = 1:length(rels)
       #   M = GAP.Globals.MappedWord(rels[i], GapObj(gens_and_invsF), GapObj(matrices_and_invs))
       #   if !isone(M)
       #      error("Group is not finite")
       #   end
       #end
        
       G2 = G.X
        
       gapMatrices = GAP.Globals.IdentityMat(length(matrices))
       for i = 1:length(matrices)
           gapMatrices[i] = Oscar.MatrixGroups._wrap_for_gap(matrices[i])
       end
       G = GAP.Globals.Group(gapMatrices)
       
       JuliaGAPMap = GAP.Globals.GroupHomomorphismByImagesNC(G,G2,GAP.Globals.GeneratorsOfGroup(G2))
       
       GAP.Globals.SetNiceMonomorphism(G,JuliaGAPMap);
       GAP.Globals.SetIsHandledByNiceMonomorphism(G, true);
       
       return G
end

end #module MatrixGroups

using .MatrixGroups
