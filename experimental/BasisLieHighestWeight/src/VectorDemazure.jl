using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted

function demazure_vw(
    L::LieAlgebraStructure, 
    reduced_expression::Vector{Int},
    lambda::Vector{ZZRingElem},
    matrices_of_operators,
)
    # Calculate vw for reduced expression w.

    # w = W([i_1, \ldots, i_s])
    # vw = f_{\alpha_{i_1}}^{l_1} \cdots f_{\alpha_{i_s}}^{l_s}.v0_

    # where l_i is defined through
    # for j in 0:s-1
    #   wj = W([i_{s+1-j},i_{s+2-j}, \ldots, i_s])
    #   l_{s-j} =  coefficients(wj*lambda)[i_{s_j}]
    # end
    R = root_system_gap(L)

    simple_roots = GAP.Globals.SimpleSystem(R)
    sparse_cartan_matrix = GAP.Globals.SparseCartanMatrix(GAPWrap.WeylGroup(R))
    s = length(reduced_expression)

    lambda_gap = GAP.Obj(Int.(lambda))

    # Calculate l_j
    l_j = []
    extremal_weight = []
    for j in 0:(s-1)
        wj_lambda = copy(lambda_gap)
        for i in s:-1:(s + 1 - j)
            GAP.Globals.ApplySimpleReflection(sparse_cartan_matrix, reduced_expression[i], wj_lambda)
        end
        
        l_smj = wj_lambda[reduced_expression[s - j]]

        pushfirst!(l_j, l_smj)

        extremal_weight = ZZ.(wj_lambda)  
    end

    # Calculate vw by applying f_{alpha_i}
    vw = sparse_row(ZZ, [(1, 1)])
    for j in length(l_j):-1:1
        for _ in 1:l_j[j]
            vw = vw * matrices_of_operators[j]
        end
    end

    return vw, extremal_weight
end
