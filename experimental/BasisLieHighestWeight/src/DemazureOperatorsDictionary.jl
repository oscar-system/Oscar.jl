function demazure_scalar_prod(beta::Int, lambda::Vector{Int})
    # only valid for simple roots
    # transpose of cartan-matrix
    # beta^v = weylgroup element w that maps beta to simple root and then let w operatore in dual space
    # coroot of a root
    return lambda[beta]
end

function demazure_s(beta::Int, lambda::Vector{Int})
    new_lambda = copy(lambda)
    new_lambda[beta] = 0
    return new_lambda
end

# Function to apply the demazure operator to a lambda vector
function demazure_operator_monom_dictionary(
    beta::Int,
    beta_wi::Vector{Int},
    lambda_old::Vector{Int}
)
    """
    Apply recursive formula
    """
    lambda = copy(lambda_old)
    scalar_prod = demazure_scalar_prod(beta, lambda)

    result = Dict{Vector{Int}, Int}()
    if scalar_prod >= 0
        for i in 0:lambda[beta]
            result[copy(lambda)] = 1
            lambda -= beta_wi
        end
    elseif scalar_prod == -1
        return result
    else
        lambda += beta_wi
        for i in 0:(-lambda[beta])
            result[copy(lambda)] = - 1
            lambda += beta_wi
        end
    end
    return result
end

# Function to update the dictionary with the demazure operator
function demazure_operator_dictionary(
    beta::Int,
    beta_wi::Vector{Int},
    dict::Dict{Vector{Int}, Int}
)
    new_dict = Dict{Vector{Int}, Int}()
    for (key, value) in dict
        updated_dict = demazure_operator_monom_dictionary(beta, beta_wi, key)
        for (new_key, new_value) in updated_dict
            # new_dict[new_key] = get(new_dict, new_key, 0) + new_value
            new_dict[new_key] = get(new_dict, new_key, 0) + value*new_value
        end
    end
    return new_dict
end


function demazure_operators_summary_dictionary(
    type::Symbol,   
    rank::Int,
    lambda::Vector{Int},
    weyl_word::Vector{Int}
)
    """
    Calculates the dimension of all weightspaces and returns them as a dictionary key => weight.
    Uses the formulas of to build up recursion.

    @misc{fourier2005tensor,
    title={Tensor product structure of affine Demazure modules and limit constructions}, 
    author={Ghislain Fourier and Peter Littelmann},
    year={2005},
    eprint={math/0412432},
    archivePrefix={arXiv},
    primaryClass={math.RT}
    }

    Start of chapter 2.2 Demazure operators(page 10)
    """
    # Initialization
    dict = Dict(lambda => 1)

    alpha_wi_list = [alpha_to_w(LieAlgebraStructure(type, rank), [i == j ? QQ(1) : QQ(0) for i in 1:rank]) for j in 1:rank]
    alpha_wi_list = [Int.(x) for x in alpha_wi_list]

    sub_word = []
    sub_dicts = [dict]
    for alpha_i in reverse(weyl_word)
        dict = demazure_operator_dictionary(alpha_i, alpha_wi_list[alpha_i], dict)

        push!(sub_word, alpha_i)
        push!(sub_dicts, copy(dict))
    end

    # Process the dictionary as needed
    return sub_dicts
end

# GAP demazure character
function get_dim_weightspace_demazure(L::LieAlgebraStructure, highest_weight::Vector{ZZRingElem}, extremal_weight::Vector{ZZRingElem}, weyl_word::Vector{Int})
    """
    Maps keys of demazure_operators_summary_dictionary to ZZ.(-key .+ highest_weight)
    """
    demazure_dict =  demazure_operators_summary_dictionary(L.lie_type, L.rank, Int.(highest_weight), weyl_word)
    demazure_dict = demazure_dict[end]
    
    result = Dict{Vector{ZZRingElem}, Int}()
    for (key, value) in demazure_dict
        # Convert the key to ZZ, apply the transformation, and assign the value
        transformed_key = ZZ.(key .- extremal_weight)
        result[transformed_key] = value
    end

    return result
end