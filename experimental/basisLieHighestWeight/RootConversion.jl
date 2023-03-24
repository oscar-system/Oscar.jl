#using Gapjm
using Oscar

fromGap = Oscar.GAP.gap_to_julia

############################################
#     conversion weight representation     #
############################################


function w_to_eps(t, n, weight)
    #println("w_to_eps: ", t, n, weight)
    if t == "A"
        return w_to_eps_A(n, weight)
    elseif t in ["B", "C", "D", "E", "F", "G"]
        #return round.(alpha_to_eps(t, n, w_to_alpha(t, n, weight))) # round not necessary
        return alpha_to_eps(t, n, w_to_alpha(t, n, weight))
    else
        println("This type of lie algebra is not supported.")
    end
end

function eps_to_w(t, n, weight)
    #println("eps_to_w: ", t, n, weight)
    if t == "A"
        return eps_to_w_A(n, weight)
    elseif t in ["B", "C", "D", "E", "F", "G"]
        return round.(alpha_to_w(t, n, eps_to_alpha(t, n, weight)))
        #return alpha_to_w(t, n, eps_to_alpha(t, n, weight))
    else
        println("This type of lie algebra is not supported.")
    end
end

function alpha_to_eps(t, n, weight)
    if t in ["B", "C", "D"]
        return alpha_to_eps_BCD(t, n, weight)
    elseif t == "E" && n in [6, 7, 8]
        return alpha_to_eps_E(n, weight)
    elseif t == "F" && n == 4
        return alpha_to_eps_F(weight)
    elseif t == "G" && n == 2
        return alpha_to_eps_G(weight)
    else
        println("This rank of lie algebra is not supported.")
    end
end

function eps_to_alpha(t, n, weight)
    if t in ["B", "C", "D"]
        return eps_to_alpha_BCD(t, n, weight)
    elseif t == "E" && n in [6, 7, 8]
        return eps_to_alpha_E(n, weight)
    elseif t == "F" && n == 4
        return eps_to_alpha_F(weight)
    elseif t == "G" && n == 2
        return eps_to_alpha_G(weight)
    else
        println("This rank of lie algebra is not supported.")
    end
end

function w_to_alpha(t, n, weight)
    C = get_CartanMatrix(t, n)
    #return round.([i for i in C*weights])
    #println("w_to_alpha: ", [i for i in C*weight])
    return [i for i in C*weight]
end

function alpha_to_w(t, n, weight)
    C_inv = get_inverse_CartanMatrix(t, n)
    #return round.([i for i in C_inv*weights])
    #println("alpha_to_w: ", [i for i in C_inv*weight])
    return [i for i in C_inv*weight]
end

function get_CartanMatrix(t, n)
    L = GAP.Globals.SimpleLieAlgebra(GAP.Obj(t), n, GAP.Globals.Rationals)
    R = GAP.Globals.RootSystem(L)
    C_list = fromGap(GAP.Globals.CartanMatrix(R))
    C = zeros(n,n)
    for i in 1:n
        for j in 1:n
            C[i,j] = C_list[i][j]
        end
    end
    return C
end

function get_inverse_CartanMatrix(t, n)
    return inv(get_CartanMatrix(t, n))
end



function alpha_to_eps_BCD(t, n, weight)
    """
    for B-D
    """
    eps = [0.0 for i=1:n]
    for i in 1:(n-1)
        eps[i] += weight[i]
        eps[i+1] -= weight[i]
    end
    if t == "B"
        eps[n] += weight[n]
    elseif t == "C"
        eps[n] += 2*weight[n]
    elseif t == "D"
        eps[n-1] += weight[n]
        eps[n] += weight[n]
    end
    return eps
end

function eps_to_alpha_BCD(t, n, weight)
    """
    for B-D
    """
    alpha = [0.0 for i=1:n]
    for i in 1:(n-2)
        alpha[i] = weight[i]
        weight[i+1] += weight[i]
    end
    if t == "B"
        alpha[n-1] = weight[n-1]
        alpha[n] += weight[n-1] + weight[n]
    elseif t == "C"
        alpha[n-1] = weight[n-1]
        alpha[n] += 0.5*weight[n-1] + 0.5*weight[n] # requires eps to be always even
    elseif t == "D"
        alpha[n-1] += (weight[n-1] - weight[n])/2
        alpha[n] += (weight[n-1] + weight[n])/2
    end
    return alpha
end

function alpha_to_eps_E(n, weight)
    """
    for E
    """
    if n == 6
        return alpha_to_eps_E6(weight)
    elseif n == 7
        return alpha_to_eps_E7(weight)
    elseif n == 8
        return alpha_to_eps_E8(weight)
    end
end

function eps_to_alpha_E(n, weight)
    """
    for E
    """
    if n == 6
        return eps_to_alpha_E6(weight)
    elseif n == 7
        return eps_to_alpha_E7(weight)
    elseif n == 8
        return eps_to_alpha_E8(weight)
    end
end

function alpha_to_eps_E6(weight)
    """
    for E6, potentially wrong order or roots (1-2-3-5-6, 3-4)
    """
    eps = [0.0 for i=1:6]
    for i=1:4
        eps[i] += weight[i]
        eps[i+1] += - weight[i]
    end
    eps[4] += weight[5]
    eps[5] += weight[5]
    for i=1:5
        eps[i] += -0.5*weight[6]
    end
    eps[6] += 0.5*sqrt(3)*weight[6]
    #println("alpha_to_eps_E6: ", eps)
    return eps
end

function eps_to_alpha_E6(weight)
    """
    for E6
    """
    alpha = [0.0 for i=1:6]
    for j=1:3
        for i=1:j
            alpha[j] += weight[i]
        end
        alpha[j] += j*(sqrt(3) / 3) *weight[6]
    end
    for i=1:4
        alpha[4] += 0.5*weight[i]
        alpha[5] += 0.5*weight[i]
    end
    alpha[4] += -0.5*weight[5] + (sqrt(3) / 2)*weight[6]
    alpha[5] += 0.5*weight[5] + 5*(sqrt(3) / 6)*weight[6]
    alpha[6] = +2*(sqrt(3) / 3)*weight[6]
    #println("eps_to_alpha_E6: ", alpha)
    return alpha
end

function alpha_to_eps_E7(weight)
    """
    for E7, potentially wrong order of roots (1-2-3-4-6-7, 4-5)
    """
    eps = [0.0 for i=1:7]
    for i=1:5
        eps[i] += weight[i]
        eps[i+1] += - weight[i]
    end
    eps[5] += weight[6]
    eps[6] += weight[6]
    for i=1:6
        eps[i] += -0.5*weight[7]
    end
    eps[7] += 0.5*sqrt(2)*weight[7]
    #println("alpha_to_eps_E7: ", eps)
    return eps
end

function eps_to_alpha_E7(weight)
    """
    for E7
    """
    alpha = [0.0 for i=1:7]
    for j=1:4
        for i=1:j
            alpha[j] += weight[i]
        end
        alpha[j] += j*(sqrt(2) / 2) *weight[7]
    end
    for i=1:5
        alpha[5] += 0.5*weight[i]
        alpha[6] += 0.5*weight[i]
    end
    alpha[5] += -0.5*weight[6] + sqrt(2)*weight[7]
    alpha[6] += 0.5*weight[6] + 3*(sqrt(2) / 2)*weight[7]
    alpha[7] = sqrt(2)*weight[7]
    #println("eps_to_alpha_E6: ", alpha)
    return alpha
end

function alpha_to_eps_E8(weight)
    """
    for E8
    """
    eps = [0.0 for i=1:8]
    for i=1:6
        eps[i] += weight[i]
        eps[i+1] += - weight[i]
    end
    eps[6] += weight[7]
    eps[7] += weight[7]
    for i=1:8
        eps[i] += -0.5*weight[8]
    end
    return eps
end

function eps_to_alpha_E8(weight)
    """
    for E8
    """
    alpha = [0.0 for i=1:8]
    for j=1:5
        for i=1:j
            alpha[j] += weight[i]
        end
        alpha[j] += -j*weight[8]
    end
    for i=1:6
        alpha[6] += 0.5*weight[i]
        alpha[7] += 0.5*weight[i]
    end
    alpha[6] += -0.5*weight[7] - 2.5*weight[8]
    alpha[7] += 0.5*weight[7] - 3.5*weight[8]
    alpha[8] = -2*weight[8]
    return alpha
end

function alpha_to_eps_F(weight) # how does this work for G?
    """
    for F
    """
    eps = [0.0 for i=1:4]
    eps[1] = weight[1] - 0.5*weight[4]
    eps[2] = - weight[1] + weight[2] - 0.5*weight[4]
    eps[3] = - weight[2] + weight[3] - 0.5*weight[4]
    eps[4] = - 0.5*weight[4]
    return eps
end

function eps_to_alpha_F(weight)
    """
    for F
    """
    alpha = [0 for i=1:4]
    alpha[1] = weight[1] - weight[4]
    alpha[2] = weight[1] + weight[2] - 2*weight[4]
    alpha[3] = weight[1] + weight[2] + weight[3] - 3*weight[4]
    alpha[4] = -2*weight[4]
    return alpha
end

function alpha_to_eps_G(weight) # how does this work for G?
    """
    for G_2
    """
    eps = [0.0 for i=1:3]
    eps[1] = weight[1] - weight[2]
    eps[2] = - weight[1] + 2*weight[2]
    eps[3] = - weight[2]
    choose_representant_eps(eps)
    return eps
end

function eps_to_alpha_G(weight)
    """
    for G_2
    """
    alpha = [0.0 for i=1:2]
    #choose_representant_eps(weight)
    if length(weight) >= 3
        weight .-= weight[3]
    end
    alpha[1] = weight[1]
    alpha[2] = (weight[1] + weight[2]) / 3
    return alpha
end

function w_to_eps_A(n, weight)
    """
    input: weight in w_i
    output: weight in eps_i
    inverse to eps_to_w
    # TODO (right now only for t=A)
    """
    res = [0 for i=1:(n+1)]
    for i in 1:n
        for l in 1:i
            res[l] += weight[i]
        end
    end
    choose_representant_eps(res)
    return res
end

function choose_representant_eps(weight)
    # choose representant eps_1 + ... + eps_m = 0
    if any(<(0), weight) # non negative
        weight .-= min(weight ...)
    end
end

function eps_to_w_A(n, weight)
    """
    input: weight in eps_i
    output: weight in w_i
    inverse to w_to_eps
    # TODO (right now only for t=A)
    """
    m = length(weight)
    choose_representant_eps(weight)
    if weight[n+1] != 0 # eps_n+1 = 0
        weight .-= weight[n+1]
        weight[n+1] = 0
    end
    #if all(>(0), weight) # minimal
    #    weight .-= min(weight ...)
    #end
    res = [0 for i=1:n]
    for i in n:-1:1
        res[i] = weight[i]
        for l in 1:i
            weight[l] -= weight[i]
        end
    end
    return res
end
