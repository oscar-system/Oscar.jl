function w_to_eps(type::String, rank::Int, weight_w::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    converts weight in rootsystem w_i to eps_i
    """
    if type in ["A", "B", "C", "D", "E", "F", "G"]
        return alpha_to_eps(type, rank, w_to_alpha(type, rank, weight_w))
    else
        println("Type needs to be one of A-D")
    end
end

function eps_to_w(type::String, rank::Int, weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    converts weight in rootsystem eps_i to w_i
    """
    if type in ["A", "B", "C", "D", "E", "F", "G"]
        # return round.(alpha_to_w(type, rank, eps_to_alpha(type, rank, weight_eps)))
        return alpha_to_w(type, rank, eps_to_alpha(type, rank, weight_eps))
    else
        println("Type needs to be one of A-D")
    end
end

function alpha_to_eps(type::String, rank::Int, weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    converts weight_alpha in rootsystem alpha_i to eps_i
    """
    if type == "A"
        return alpha_to_eps_A(rank, weight_alpha)
    elseif type in ["B", "C", "D"]
        return alpha_to_eps_BCD(type, rank, weight_alpha)
    elseif type == "E" && rank in [6, 7, 8]
        return alpha_to_eps_E(rank, weight_alpha)
    elseif type == "F" && rank == 4
        return alpha_to_eps_F(weight_alpha)
    elseif type == "G" && rank == 2
        return alpha_to_eps_G(weight_alpha)
    else
        println("This rank of lie algebra is not supported.")
    end
end

function eps_to_alpha(type::String, rank::Int, weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    converts weight_eps in rootsystem eps_i to alpha_i
    """
    if type == "A"
        return eps_to_alpha_A(rank, weight_eps)
    elseif type in ["B", "C", "D"]
        return eps_to_alpha_BCD(type, rank, weight_eps)
    elseif type == "E" && rank in [6, 7, 8]
        return eps_to_alpha_E(rank, weight_eps)
    elseif type == "F" && rank == 4
        return eps_to_alpha_F(weight_eps)
    elseif type == "G" && rank == 2
        return eps_to_alpha_G(weight_eps)
    else
        println("This rank of lie algebra is not supported.")
    end
end

function w_to_alpha(type, rank, weight_w::Vector{QQFieldElem})::Vector{QQFieldElem}
    C = get_inverse_CartanMatrix(type, rank)
    return [i for i in C*weight_w]
end

function alpha_to_w(type::String, rank::Int, weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    C_inv = get_CartanMatrix(type, rank)
    return [i for i in C_inv*weight_alpha]
end

function get_CartanMatrix(type::String, rank::Int)::Matrix{QQFieldElem}
    L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    R = GAP.Globals.RootSystem(L)
    C = Matrix{QQFieldElem}(GAP.Globals.CartanMatrix(R))
    # println("C: ", C)
    return C
end

function get_inverse_CartanMatrix(type::String, rank::Int)::Matrix{QQFieldElem}
    return inv(get_CartanMatrix(type, rank))
end

function alpha_to_eps_BCD(type::String, rank::Int, weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for B-D
    """
    weight_eps = [QQ(0) for i in 1:rank]
    for i in 1:(rank - 1)
        weight_eps[i] += weight_alpha[i]
        weight_eps[i+1] -= weight_alpha[i]
    end
    if type == "B"
        weight_eps[rank] += weight_alpha[rank]
    elseif type == "C"
        weight_eps[rank] += QQ(2)*weight_alpha[rank]
    elseif type == "D"
        weight_eps[rank - 1] += weight_alpha[rank]
        weight_eps[rank] += weight_alpha[rank]
    end
    return weight_eps
end

function eps_to_alpha_BCD(type::String, rank::Int, weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for B-D
    """
    weight_alpha = [QQ(0) for i in 1:rank]
    for i in 1:(rank - 2)
        weight_alpha[i] = weight_eps[i]
        weight_eps[i+1] += weight_eps[i]
    end
    if type == "B"
        weight_alpha[rank - 1] = weight_eps[rank - 1]
        weight_alpha[rank] += weight_eps[rank-1] + weight_eps[rank]
    elseif type == "C"
        weight_alpha[rank - 1] = weight_eps[rank - 1]
        weight_alpha[rank] += QQ(1, 2)*weight_eps[rank - 1] + QQ(1, 2)*weight_eps[rank]
    elseif type == "D"
        weight_alpha[rank - 1] += (weight_eps[rank - 1] - QQ(1, 2)*weight_eps[rank])
        weight_alpha[rank] += (weight_eps[rank - 1] + QQ(1, 2)*weight_eps[rank])
    end
    return weight_alpha
end

function alpha_to_eps_E(rank::Int, weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for E
    """
    if rank == 6
        return alpha_to_eps_E6(weight_alpha)
    elseif rank == 7
        return alpha_to_eps_E7(weight_alpha)
    elseif rank == 8
        return alpha_to_eps_E8(weight_alpha)
    end
end

function eps_to_alpha_E(rank::Int, weight_eps)
    """
    for E
    """
    if rank == 6
        return eps_to_alpha_E6(weight_eps)
    elseif rank == 7
        return eps_to_alpha_E7(weight_eps)
    elseif rank == 8
        return eps_to_alpha_E8(weight_eps)
    end
end

function alpha_to_eps_E6(weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for E6, potentially wrong order or roots (1-2-3-5-6, 3-4)
    """
    weight_eps = [QQ(0) for i in 1:6]
    for i in 1:4
        weight_eps[i] += weight_alpha[i]
        weight_eps[i + 1] += - weight_alpha[i]
    end
    weight_eps[4] += weight_alpha[5]
    weight_eps[5] += weight_alpha[5]
    for i in 1:5
        weight_eps[i] += QQ(-1, 2)*weight_alpha[6]
    end
    weight_eps[6] += QQ(1, 2)*sqrt(3)*weight_alpha[6]
    return eps
end

function eps_to_alpha_E6(weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem} # TODO sqrt
    """
    for E6
    """
    weight_alpha = [QQ(0) for i in 1:6]
    for j in 1:3
        for i in 1:j
            weight_alpha[j] += weight_eps[i]
        end
        weight_alpha[j] += j*(sqrt(3) / 3) *weight_eps[6]
    end
    for i in 1:4
        weight_alpha[4] += 0.5*weight_eps[i]
        weight_alpha[5] += 0.5*weight_eps[i]
    end
    weight_alpha[4] += -0.5*weight_eps[5] + (sqrt(3) / 2)*weight_eps[6]
    weight_alpha[5] += 0.5*weight_eps[5] + 5*(sqrt(3) / 6)*weight_eps[6]
    weight_alpha[6] = +2*(sqrt(3) / 3)*weight_eps[6]
    #println("eps_to_alpha_E6: ", alpha)
    return weight_alpha
end

function alpha_to_eps_E7(weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem} # TODO sqrt
    """
    for E7, potentially wrong order of roots (1-2-3-4-6-7, 4-5)
    """
    weight_eps = [QQ(0) for i in 1:7]
    for i in 1:5
        weight_eps[i] += weight_alpha[i]
        weight_eps[i + 1] += - weight_alpha[i]
    end
    weight_eps[5] += weight_alpha[6]
    weight_eps[6] += weight_alpha[6]
    for i in 1:6
        weight_eps[i] += QQ(-1, 2)*weight_alpha[7]
    end
    weight_eps[7] += 0.5*sqrt(2)*weight_alpha[7]
    return weight_eps
end

function eps_to_alpha_E7(weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem} # TODO sqrt
    """
    for E7
    """
    weight_alpha = [QQ(0) for i in 1:7]
    for j in 1:4
        for i in 1:j
            weight_alpha[j] += weight_eps[i]
        end
        weight_alpha[j] += j*(sqrt(2) / 2) *weight_eps[7]
    end
    for i in 1:5
        weight_alpha[5] += 0.5*weight_eps[i]
        weight_alpha[6] += 0.5*weight_eps[i]
    end
    weight_alpha[5] += -0.5*weight_eps[6] + sqrt(2)*weight_eps[7]
    weight_alpha[6] += 0.5*weight_eps[6] + 3*(sqrt(2) / 2)*weight_eps[7]
    weight_alpha[7] = sqrt(2)*weight_eps[7]
    return weight_alpha
end

function alpha_to_eps_E8(weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for E8
    """
    weight_eps = [0.0 for i in 1:8]
    for i in 1:6
        weight_eps[i] += weight_alpha[i]
        weight_eps[i+1] += - weight_alpha[i]
    end
    weight_eps[6] += weight_alpha[7]
    weight_eps[7] += weight_alpha[7]
    for i in 1:8
        weight_eps[i] += -0.5*weight_alpha[8]
    end
    return weight_eps
end

function eps_to_alpha_E8(weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for E8
    """
    weight_alpha = [QQ(0) for i in 1:8]
    for j in 1:5
        for i in 1:j
            weight_alpha[j] += weight_eps[i]
        end
        weight_alpha[j] += -j*weight_eps[8]
    end
    for i in 1:6
        weight_alpha[6] += QQ(1, 2)*weight_eps[i]
        weight_alpha[7] += QQ(1, 2)*weight_eps[i]
    end
    weight_alpha[6] += QQ(-1, 2)*weight_eps[7] - QQ(5, 2)*weight_eps[8]
    weight_alpha[7] += QQ(1, 2)*weight_eps[7] - QQ(7, 2)*weight_eps[8]
    weight_alpha[8] = QQ(-2)*weight_eps[8]
    return alpha
end

function alpha_to_eps_F(weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for F
    """
    weight_eps = [QQ(0) for i in 1:4]
    weight_eps[1] = weight_alpha[1] - QQ(1, 2)*weight_alpha[4]
    weight_eps[2] = - weight_alpha[1] + weight_alpha[2] - QQ(1, 2)*weight_alpha[4]
    weight_eps[3] = - weight_alpha[2] + weight_alpha[3] - QQ(1, 2)*weight_alpha[4]
    weight_eps[4] = - QQ(1, 2)*weight_alpha[4]
    return weight_eps
end

function eps_to_alpha_F(weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for F
    """
    weight_alpha = [QQ(0) for i in 1:4]
    weight_alpha[1] = weight_eps[1] - weight_eps[4]
    weight_alpha[2] = weight_eps[1] + weight_eps[2] - QQ(2)*weight_eps[4]
    weight_alpha[3] = weight_eps[1] + weight_eps[2] + weight_eps[3] - QQ(3)*weight_eps[4]
    weight_alpha[4] = QQ(-2)*weight_eps[4]
    return weight_alpha
end

function alpha_to_eps_G(weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for G_2
    """
    weight_eps = [QQ(0) for i in 1:3]
    weight_eps[1] = weight_alpha[1] - weight_alpha[2]
    weight_eps[2] = - weight_alpha[1] + QQ(2)*weight_alpha[2]
    weight_eps[3] = - weight_alpha[2]
    choose_representant_eps(weight_eps)
    return weight_eps
end

function eps_to_alpha_G(weight_eps::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for G_2
    """
    weight_alpha = [QQ(0) for i in 1:2]
    if length(weight_eps) >= 3
        weight_eps .-= weight_eps[3]
    end
    weight_alpha[1] = weight_eps[1]
    weight_alpha[2] = QQ(1, 3)*(weight_eps[1] + weight_eps[2])
    return weight_alpha
end

function choose_representant_eps(weight_eps::Vector{QQFieldElem})
    # choose representant eps_1 + ... + eps_m = 0
    if any(<(0), weight_eps) # non negative
        weight_eps .-= min(weight_eps ...)
    end
end

function alpha_to_eps_A(rank::Int, weight_alpha::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for A
    """
    weight_eps = [QQ(0) for i in 1:(rank + 1)]
    for i in 1:rank
        weight_eps[i] += weight_alpha[i]
        weight_eps[i + 1] -= weight_alpha[i]
    end
    choose_representant_eps(weight_eps)
    return weight_eps
end

function eps_to_alpha_A(rank::Int, weight_eps_original::Vector{QQFieldElem})::Vector{QQFieldElem}
    """
    for A
    """
    weight_eps = copy(weight_eps_original) # Original weight should not get modified
    if length(weight_eps) == rank
        push!(weight_eps, QQ(0))
    end
    weight_alpha = [QQ(0) for i in 1:(rank + 1)]
    for i in 1:(rank + 1)
        for j in 1:i
            weight_alpha[i] += weight_eps[j]
        end
    end
    m = weight_alpha[rank + 1]*QQ(1, rank + 1)
    for i in 1:rank
        weight_alpha[i] -= QQ(i)*m
    end
    pop!(weight_alpha)
    return weight_alpha
end
