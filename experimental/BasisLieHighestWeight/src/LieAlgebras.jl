# ?

using Oscar
#using SparseArrays

fromGap = Oscar.GAP.gap_to_julia


function lieAlgebra(t::String, n::Int)
    """
    Creates the Lie-algebra as a GAP object that gets used for a lot other computations with GAP
    """
    L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(t), n, GAP.Globals.Rationals)
    return L, GAP.Globals.ChevalleyBasis(L)
end


gapReshape(A) = sparse_matrix(QQ, hcat(A...))
#gapReshape(A) = sparse(hcat(A...))

# temporary workaround for issue 2128
function multiply_scalar(A::SMat{T}, d) where T
    for i in 1:nrows(A)
        scale_row!(A, i, T(d))
    end
    return A
    #return identity_matrix(SMat, QQ, size(A)[1])*A
end

function matricesForOperators(L, hw, ops)
    """
    used to create tensorMatricesForOperators
    """
    #M = GAP.Globals.HighestWeightModule(L, GAP.Obj(hw))
    #mats = GAP.Globals.List(ops, o -> GAP.Globals.MatrixOfAction(GAP.Globals.Basis(M), o))
    M = Oscar.GAP.Globals.HighestWeightModule(L, Oscar.GAP.julia_to_gap(hw))
    mats = Oscar.GAP.Globals.List(ops, o -> Oscar.GAP.Globals.MatrixOfAction(GAP.Globals.Basis(M), o))
    #mats = gapReshape.(fromGap(mats))
    mats = gapReshape.( Oscar.GAP.gap_to_julia(mats))
    denominators = map(y->denominator(y[2]), union(union(mats...)...))
    #d = convert(QQ, lcm(denominators))
    d = lcm(denominators)# // 1
    mats = (A->change_base_ring(ZZ, multiply_scalar(A, d))).(mats)
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
