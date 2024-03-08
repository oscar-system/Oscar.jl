# This file implements explicit models of symplectic reflection groups.
#
# Ulrich Thiel, 2024


function symplectic_reflection_group(W::MatrixGroup)

    if !is_complex_reflection_group(W)
        throw(ArgumentError("Group must be a complex reflection group"))
    end

    K = base_ring(W)
    n = degree(W)
    matspace = matrix_space(K, 2*n, 2*n)
    gens_symp = elem_type(matspace)[]

    for w in gens(W)
        w_symp = matspace(block_diagonal_matrix([matrix(w), transpose(matrix(w^-1))]))
        push!(gens_symp, w_symp)
    end

    W_symp = matrix_group(gens_symp)

    set_attribute!(W_symp, :order, order(W))
    set_attribute!(W_symp, :is_symplectic_reflection_group, true)

    return W_symp
end

function is_symplectic_reflection_group(G::MatrixGroup)
    if has_attribute(G, :is_symplectic_reflection_group)
        return get_attribute(G, :is_symplectic_reflection_group)
    end
    return false 
    #this should be upgraded later to work with a general matrix group
end
