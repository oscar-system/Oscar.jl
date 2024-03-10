#R = root_system(:X, n)
#W = weyl_group(R)

#Sie brauchen eine reduzierte Zerlegung von w. Vielleicht können Sie die zu Beginn als Input annehmen
#w = W([i_1, \ldots, i_s])
#Dazu haben Sie lambda = WeightLatticeElem(R,[m_1, \ldots, m_n]) als höchstes Gewicht.

#Ihr Vektor ist
#vw = f_{\alpha_{i_1}}^{l_1} \cdots f_{\alpha_{i_s}}^{l_s}.v0_

#hierbei bestimmen wir l_i durch
#for j in 0:s-1
#  wj = W([i_{s+1-j},i_{s+2-j}, \ldots, i_s])
#  l_{s-j} =  coefficients(wj*lambda)[i_{s_j}]  # TODO >= 0
#end


using ..Oscar
using ..Oscar: GAPWrap, IntegerUnion, _is_weighted

function demazure_vw(
    L::LieAlgebraStructure, 
    reduced_expression::Vector{Int},
    lambda::Vector{ZZRingElem},
    matrices_of_operators,
)
    R = root_system_gap(L)

    simple_roots = GAP.Globals.SimpleSystem(R)
    sparse_cartan_matrix = GAP.Globals.SparseCartanMatrix(GAPWrap.WeylGroup(R))
    s = length(reduced_expression)

    # Write lambda as sum of simple roots in GAP
    lambda_gap = simple_roots[1] - simple_roots[1]
    for i in 1:length(lambda)
        lambda_gap = lambda_gap + Int(lambda[i]) * simple_roots[i]
    end
    println("lambda in GAP: ", lambda_gap)

    # Calculate l_j = w^{-1} w_{s+1-j} ... w^{-1} w_{s} lambda
    l_j = []
    for j in 0:(s-1)
        println("")
        println("j = ", j)
        wj_lambda = copy(lambda_gap)
        for i in s:-1:(s + 1 - j)
            println("i = ", i)
            GAP.Globals.ApplySimpleReflection(sparse_cartan_matrix, reduced_expression[i], wj_lambda)
        end

        println("start twist: ", wj_lambda)
        # Twist with w^{-1} = w_s ... w_1
        for k in 1:s
            print(" ", k)
            GAP.Globals.ApplySimpleReflection(sparse_cartan_matrix, reduced_expression[k], wj_lambda)
        end
        #println("middle twist: ", wj_lambda)
        #for k in 1:s
        #    print(" ", s-k+1)
        #    GAP.Globals.ApplySimpleReflection(sparse_cartan_matrix, reduced_expression[s-k+1], wj_lambda)
        #end
        println("end twist: ", wj_lambda)

        l_smj = Int.(wj_lambda)[reduced_expression[s - j]]
        println("w_j * lambda = ", wj_lambda)
        println("l_{s-j} = ", l_smj)
        pushfirst!(l_j, l_smj)
    end
    println("l_j: ", l_j)

    # Calculate vw
    vw = sparse_row(ZZ, [(1, 1)])
    for j in length(l_j):-1:1
        println("j: ", j, " l_j:", l_j[j])
        for _ in 1:l_j[j]
            vw = vw * matrices_of_operators[j]
        end
        println("vw, j: ", j, " ", vw)
    end

    println("vw: ", vw)

    return vw
end

#lie_alg = lie_algebra(:A, 2)
#reduced_expression = Vector{Int}([1, 2, 1])
#lambda = [1, 1]
#matrices_of_operators = tensor_matrices_of_operators(
#    L, ZZ.(lambda), birational_sequence.operators
#  )

#println(demazure_vw(lie_alg, reduced_expression, lambda, matrices_of_operators))
