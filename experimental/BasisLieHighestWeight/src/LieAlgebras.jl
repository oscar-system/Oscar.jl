fromGap = Oscar.GAP.gap_to_julia


function create_lie_algebra(type::String, rank::Int)::GAP.Obj
    """
    Creates the Lie-algebra as a GAP object that gets used for a lot of other computations with GAP
    """
    lie_algebra = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    return lie_algebra
end


gapReshape(A) = sparse_matrix(QQ, hcat(A...))

# temporary workaround for issue 2128
function multiply_scalar(A::SMat{T}, d) where T
    for i in 1:nrows(A)
        scale_row!(A, i, T(d))
    end
    return A
end

function matricesForOperators(lie_algebra::GAP.Obj, highest_weight::Vector{ZZRingElem}, 
                              operators::GAP.Obj)::Vector{SMat{ZZRingElem}}
    """
    used to create tensorMatricesForOperators
    """
    highest_weight_int = convert(Vector{Int}, highest_weight)
    M = GAP.Globals.HighestWeightModule(lie_algebra, GAP.julia_to_gap(highest_weight_int))
    matrices_of_operators = GAP.Globals.List(operators, o -> GAP.Globals.MatrixOfAction(GAPWrap.Basis(M), o))
    matrices_of_operators = gapReshape.(GAP.gap_to_julia(matrices_of_operators))
    denominators = map(y->denominator(y[2]), union(union(matrices_of_operators...)...))
    common_denominator = lcm(denominators)# // 1
    matrices_of_operators = (A->change_base_ring(ZZ, multiply_scalar(A, common_denominator))).(matrices_of_operators)
    return matrices_of_operators
end


function weights_for_operators(lie_algebra::GAP.Obj, cartan::GAP.Obj, operators::GAP.Obj)::Vector{Vector{ZZRingElem}}
    """
    Calculates the weight weights[i] in w_i for each operator operators[i]
    """
    """cartan = [Vector{Int}(x) for x in GAP.Globals.ExtRepOfObj.(cartan)]
    operators = [Vector{Int}(x) for x in GAP.Globals.ExtRepOfObj.(operators)]#
    if any(iszero.(operators))
        error("ops should be non-zero")
    end
    println([findfirst(v .!= 0) for v in operators])
    
    return [
        [(dot(h, v))[findfirst(v .!= 0)] / (v)[findfirst(v .!= 0)] for h in cartan] for v in operators
    ]


    """
    # TODO delete fromGap. Multiplication of cartan and operators is not regular matrix multiplication
    cartan = fromGap(cartan, recursive=false)
    operators = fromGap(operators, recursive=false)
    asVec(v) = fromGap(GAPWrap.ExtRepOfObj(v))
    #println(cartan)
    #println(operators)
    if any(iszero.(asVec.(operators)))
        error("ops should be non-zero")
    end
    nzi(v) = findfirst(asVec(v) .!= 0)
    #println([nzi(v) for v in operators])
    for h in cartan
        for v in operators
            #println("-")
            #println(asVec(v))
            #println(asVec(h))
            #println(asVec(h*v))
        end
    end

    return [
        [ZZ(QQ(asVec(h*v)[nzi(v)], asVec(v)[nzi(v)])) for h in cartan] for v in operators
    ]
end
