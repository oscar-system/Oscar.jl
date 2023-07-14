function w_to_eps(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    """
    converts weight in rootsystem w_i to eps_i
    """
    if type in ["A", "B", "C", "D", "E", "F", "G"]
        return alpha_to_eps(type, rank, w_to_alpha(type, rank, weight))
    else
        println("Type needs to be one of A-D")
    end
end

function eps_to_w(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    """
    converts weight in rootsystem eps_i to w_i
    """
    if type in ["A", "B", "C", "D", "E", "F", "G"]
        return round.(alpha_to_w(type, rank, eps_to_alpha(type, rank, weight)))
    else
        println("Type needs to be one of A-D")
    end
end

function alpha_to_eps(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    """
    converts weight in rootsystem alpha_i to eps_i
    """
    if type == "A"
        return alpha_to_eps_A(rank, weight)
    elseif type in ["B", "C", "D"]
        return alpha_to_eps_BCD(type, rank, weight)
    elseif type == "E" && rank in [6, 7, 8]
        return alpha_to_eps_E(rank, weight)
    elseif type == "F" && rank == 4
        return alpha_to_eps_F(weight)
    elseif type == "G" && rank == 2
        return alpha_to_eps_G(weight)
    else
        println("This rank of lie algebra is not supported.")
    end
end

function eps_to_alpha(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    """
    converts weight in rootsystem eps_i to alpha_i
    """
    if type == "A"
        return eps_to_alpha_A(rank, weight)
    elseif type in ["B", "C", "D"]
        return eps_to_alpha_BCD(type, rank, weight)
    elseif type == "E" && rank in [6, 7, 8]
        return eps_to_alpha_E(rank, weight)
    elseif type == "F" && rank == 4
        return eps_to_alpha_F(weight)
    elseif type == "G" && rank == 2
        return eps_to_alpha_G(weight)
    else
        println("This rank of lie algebra is not supported.")
    end
end

function w_to_alpha(type, rank, weight::Vector{Int})::Vector{Int}
    C = get_CartanMatrix(type, rank)
    return [i for i in C*weight]
end

function alpha_to_w(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    C_inv = get_inverse_CartanMatrix(type, rank)
    return [i for i in C_inv*weight]
end

function get_CartanMatrix(type::String, rank::Int)
    L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(type), rank, GAP.Globals.Rationals)
    R = GAP.Globals.RootSystem(L)
    C = Matrix{Int}(GAP.Globals.CartanMatrix(R))
    return C
end

function get_inverse_CartanMatrix(type::String, rank::Int)
    return inv(get_CartanMatrix(type, rank))
end

function alpha_to_eps_BCD(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    """
    for B-D
    """
    eps = [0.0 for i in 1:rank]
    for i in 1:(rank-1)
        eps[i] += weight[i]
        eps[i+1] -= weight[i]
    end
    if type == "B"
        eps[rank] += weight[rank]
    elseif type == "C"
        eps[rank] += 2*weight[rank]
    elseif type == "D"
        eps[rank - 1] += weight[rank]
        eps[rank] += weight[rank]
    end
    return eps
end

function eps_to_alpha_BCD(type::String, rank::Int, weight::Vector{Int})::Vector{Int}
    """
    for B-D
    """
    alpha = [0.0 for i in 1:rank]
    for i in 1:(rank-2)
        alpha[i] = weight[i]
        weight[i+1] += weight[i]
    end
    if type == "B"
        alpha[rank - 1] = weight[rank - 1]
        alpha[rank] += weight[rank-1] + weight[rank]
    elseif type == "C"
        alpha[rank - 1] = weight[rank - 1]
        alpha[rank] += 0.5*weight[rank - 1] + 0.5*weight[rank] # requires eps to be always even
    elseif type == "D"
        alpha[rank - 1] += (weight[rank - 1] - weight[rank])/2
        alpha[rank] += (weight[rank - 1] + weight[rank])/2
    end
    return alpha
end

function alpha_to_eps_E(rank::Int, weight::Vector{Int})::Vector{Int}
    """
    for E
    """
    if rank == 6
        return alpha_to_eps_E6(weight)
    elseif rank == 7
        return alpha_to_eps_E7(weight)
    elseif rank == 8
        return alpha_to_eps_E8(weight)
    end
end

function eps_to_alpha_E(rank::Int, weight)
    """
    for E
    """
    if rank == 6
        return eps_to_alpha_E6(weight)
    elseif rank == 7
        return eps_to_alpha_E7(weight)
    elseif rank == 8
        return eps_to_alpha_E8(weight)
    end
end

function alpha_to_eps_E6(weight::Vector{Int})::Vector{Int}
    """
    for E6, potentially wrong order or roots (1-2-3-5-6, 3-4)
    """
    eps = [0.0 for i in 1:6]
    for i in 1:4
        eps[i] += weight[i]
        eps[i + 1] += - weight[i]
    end
    eps[4] += weight[5]
    eps[5] += weight[5]
    for i in 1:5
        eps[i] += -0.5*weight[6]
    end
    eps[6] += 0.5*sqrt(3)*weight[6]
    return eps
end

function eps_to_alpha_E6(weight)
    """
    for E6
    """
    alpha = [0.0 for i in 1:6]
    for j in 1:3
        for i in 1:j
            alpha[j] += weight[i]
        end
        alpha[j] += j*(sqrt(3) / 3) *weight[6]
    end
    for i in 1:4
        alpha[4] += 0.5*weight[i]
        alpha[5] += 0.5*weight[i]
    end
    alpha[4] += -0.5*weight[5] + (sqrt(3) / 2)*weight[6]
    alpha[5] += 0.5*weight[5] + 5*(sqrt(3) / 6)*weight[6]
    alpha[6] = +2*(sqrt(3) / 3)*weight[6]
    #println("eps_to_alpha_E6: ", alpha)
    return alpha
end

function alpha_to_eps_E7(weight::Vector{Int})::Vector{Int}
    """
    for E7, potentially wrong order of roots (1-2-3-4-6-7, 4-5)
    """
    eps = [0.0 for i in 1:7]
    for i in 1:5
        eps[i] += weight[i]
        eps[i + 1] += - weight[i]
    end
    eps[5] += weight[6]
    eps[6] += weight[6]
    for i in 1:6
        eps[i] += -0.5*weight[7]
    end
    eps[7] += 0.5*sqrt(2)*weight[7]
    return eps
end

function eps_to_alpha_E7(weight::Vector{Int})::Vector{Int}
    """
    for E7
    """
    alpha = [0.0 for i in 1:7]
    for j in 1:4
        for i in 1:j
            alpha[j] += weight[i]
        end
        alpha[j] += j*(sqrt(2) / 2) *weight[7]
    end
    for i in 1:5
        alpha[5] += 0.5*weight[i]
        alpha[6] += 0.5*weight[i]
    end
    alpha[5] += -0.5*weight[6] + sqrt(2)*weight[7]
    alpha[6] += 0.5*weight[6] + 3*(sqrt(2) / 2)*weight[7]
    alpha[7] = sqrt(2)*weight[7]
    return alpha
end

function alpha_to_eps_E8(weight::Vector{Int})::Vector{Int}
    """
    for E8
    """
    eps = [0.0 for i in 1:8]
    for i in 1:6
        eps[i] += weight[i]
        eps[i+1] += - weight[i]
    end
    eps[6] += weight[7]
    eps[7] += weight[7]
    for i in 1:8
        eps[i] += -0.5*weight[8]
    end
    return eps
end

function eps_to_alpha_E8(weight::Vector{Int})::Vector{Int}
    """
    for E8
    """
    alpha = [0.0 for i in 1:8]
    for j in 1:5
        for i in 1:j
            alpha[j] += weight[i]
        end
        alpha[j] += -j*weight[8]
    end
    for i in 1:6
        alpha[6] += 0.5*weight[i]
        alpha[7] += 0.5*weight[i]
    end
    alpha[6] += -0.5*weight[7] - 2.5*weight[8]
    alpha[7] += 0.5*weight[7] - 3.5*weight[8]
    alpha[8] = -2*weight[8]
    return alpha
end

function alpha_to_eps_F(weight::Vector{Int})::Vector{Int}
    """
    for F
    """
    eps = [0.0 for i in 1:4]
    eps[1] = weight[1] - 0.5*weight[4]
    eps[2] = - weight[1] + weight[2] - 0.5*weight[4]
    eps[3] = - weight[2] + weight[3] - 0.5*weight[4]
    eps[4] = - 0.5*weight[4]
    return eps
end

function eps_to_alpha_F(weight::Vector{Int})::Vector{Int}
    """
    for F
    """
    alpha = [0 for i in 1:4]
    alpha[1] = weight[1] - weight[4]
    alpha[2] = weight[1] + weight[2] - 2*weight[4]
    alpha[3] = weight[1] + weight[2] + weight[3] - 3*weight[4]
    alpha[4] = -2*weight[4]
    return alpha
end

function alpha_to_eps_G(weight::Vector{Int})::Vector{Int}
    """
    for G_2
    """
    eps = [0.0 for i in 1:3]
    eps[1] = weight[1] - weight[2]
    eps[2] = - weight[1] + 2*weight[2]
    eps[3] = - weight[2]
    choose_representant_eps(eps)
    return eps
end

function eps_to_alpha_G(weight::Vector{Int})::Vector{Int}
    """
    for G_2
    """
    alpha = [0.0 for i in 1:2]
    if length(weight) >= 3
        weight .-= weight[3]
    end
    alpha[1] = weight[1]
    alpha[2] = (weight[1] + weight[2]) / 3
    return alpha
end

function choose_representant_eps(weight::Vector{Int})
    # choose representant eps_1 + ... + eps_m = 0
    if any(<(0), weight) # non negative
        weight .-= min(weight ...)
    end
end

function alpha_to_eps_A(rank::Int, weight::Vector{Int})::Vector{Int}
    """
    for A
    """
    eps = [0 for i in 1:(rank + 1)]
    for i in 1:rank
        eps[i] += weight[i]
        eps[i + 1] -= weight[i]
    end
    choose_representant_eps(eps)
    return eps
end

function eps_to_alpha_A(rank::Int, weight::Vector{Int})::Vector{Int}
    """
    for A
    """
    if length(weight) == rank
        append!(weight, 0)
    end
    alpha = [0.0 for i in 1:(rank + 1)]
    for i in 1:(rank + 1)
        for j in 1:i
            alpha[i] += weight[j]
        end
    end
    m = alpha[rank + 1] / (rank + 1)
    for i in 1:rank
        alpha[i] -= i*m
    end
    pop!(alpha)
    return alpha
end
