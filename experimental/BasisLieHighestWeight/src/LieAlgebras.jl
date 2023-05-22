using Oscar

fromGap = Oscar.GAP.gap_to_julia


function create_lie_lgebra(type::String, rank::Int)::Tuple{GAP.Obj, GAP.Obj}
    """
    Creates the Lie-algebra as a GAP object that gets used for a lot of other computations with GAP
    """
    lie_algebra = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    return lie_algebra, GAP.Globals.ChevalleyBasis(lie_algebra)
end


gapReshape(A) = sparse_matrix(QQ, hcat(A...))

# temporary workaround for issue 2128
function multiply_scalar(A::SMat{T}, d) where T
    for i in 1:nrows(A)
        scale_row!(A, i, T(d))
    end
    return A
end

function matricesForOperators(lie_algebra::GAP.Obj, highest_weight::Vector{Int}, 
                              ops::GAP.Obj)::Vector{SMat{ZZRingElem}}
    """
    used to create tensorMatricesForOperators
    """
    M = Oscar.GAP.Globals.HighestWeightModule(lie_algebra, Oscar.GAP.julia_to_gap(highest_weight))
    matrices_of_operators = Oscar.GAP.Globals.List(ops, o -> Oscar.GAP.Globals.MatrixOfAction(GAP.Globals.Basis(M), o))
    matrices_of_operators = gapReshape.( Oscar.GAP.gap_to_julia(matrices_of_operators))
    denominators = map(y->denominator(y[2]), union(union(matrices_of_operators...)...))
    common_denominator = lcm(denominators)# // 1
    matrices_of_operators = (A->change_base_ring(ZZ, multiply_scalar(A, common_denominator))).(matrices_of_operators)
    return matrices_of_operators
end


function weights_for_operators(lie_algebra::GAP.Obj, cartan::GAP.Obj, operators::GAP.Obj)::Vector{Vector{Int}}
    """
    Calculates the weight wts[i] for each operator ops[i]
    """
    cartan = fromGap(cartan, recursive=false)
    operators = fromGap(operators, recursive=false)
    asVec(v) = fromGap(GAP.Globals.ExtRepOfObj(v))
    if any(iszero.(asVec.(operators)))
        error("ops should be non-zero")
    end
    nzi(v) = findfirst(asVec(v) .!= 0)
    return [
        [asVec(h*v)[nzi(v)] / asVec(v)[nzi(v)] for h in cartan] for v in operators
    ]
end
