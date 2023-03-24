# ?

using Oscar
using SparseArrays

G = Oscar.GAP.Globals
forGap = Oscar.GAP.julia_to_gap
fromGap = Oscar.GAP.gap_to_julia


function lieAlgebra(t::String, n::Int)
    """
    Creates the Lie-algebra as a GAP object that gets used for a lot other computations with GAP
    """
    L = G.SimpleLieAlgebra(forGap(t), n, G.Rationals)
    return L, G.ChevalleyBasis(L)
end


gapReshape(A) = sparse(hcat(A...))


function matricesForOperators(L, hw, ops)
    """
    used to create tensorMatricesForOperators
    """
    M = G.HighestWeightModule(L, forGap(hw))
    mats = G.List(ops, o -> G.MatrixOfAction(G.Basis(M), o))
    mats = gapReshape.(fromGap(mats))
    d = lcm(denominator.(union(mats...)))
    mats = (A->ZZ.(A*d)).(mats)
    return mats
end


function weightsForOperators(L, cartan, ops)
    """
    Calculates the weight wts[i] for each operator ops[i]
    """
    cartan = fromGap(cartan, recursive=false)
    ops = fromGap(ops, recursive=false)
    asVec(v) = fromGap(G.ExtRepOfObj(v))
    if any(iszero.(asVec.(ops)))
        error("ops should be non-zero")
    end
    nzi(v) = findfirst(asVec(v) .!= 0)
    return [
        [asVec(h*v)[nzi(v)] / asVec(v)[nzi(v)] for h in cartan] for v in ops
    ]
end