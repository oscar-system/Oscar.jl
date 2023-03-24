# ?

using Oscar
using SparseArrays

fromGap = Oscar.GAP.gap_to_julia


function lieAlgebra(t::String, n::Int)
    """
    Creates the Lie-algebra as a GAP object that gets used for a lot other computations with GAP
    """
    L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(t), n, GAP.Globals.Rationals)
    return L, GAP.Globals.ChevalleyBasis(L)
end


gapReshape(A) = sparse(hcat(A...))


function matricesForOperators(L, hw, ops)
    """
    used to create tensorMatricesForOperators
    """
    M = GAP.Globals.HighestWeightModule(L, GAP.Obj(hw))
    mats = GAP.Globals.List(ops, o -> GAP.Globals.MatrixOfAction(GAP.Globals.Basis(M), o))
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
    asVec(v) = fromGap(GAP.Globals.ExtRepOfObj(v))
    if any(iszero.(asVec.(ops)))
        error("ops should be non-zero")
    end
    nzi(v) = findfirst(asVec(v) .!= 0)
    return [
        [asVec(h*v)[nzi(v)] / asVec(v)[nzi(v)] for h in cartan] for v in ops
    ]
end
