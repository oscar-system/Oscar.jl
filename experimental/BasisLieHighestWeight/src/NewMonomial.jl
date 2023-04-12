include("./VectorSpaceBases.jl")
include("./TensorModels.jl")
include("./LieAlgebras.jl")
include("./MonomialOrder.jl")
include("./WeylPolytope.jl")

fromGap = Oscar.GAP.gap_to_julia


function calc_weight(mon::ZZMPolyRingElem, weights::Vector{Vector{Int}})::Vector{Int}
    """
    calculates weight associated with monomial mon
    """
    degree_mon = degrees(mon)
    weight = [0 for i in 1:length(weights[1])]
    for i in 1:length(degree_mon)
        weight .+= degree_mon[i] * weights[i]
    end
    return weight
end

function calc_vec(v0::SRow{ZZRingElem}, mon::ZZMPolyRingElem, 
                  matrices_of_operators::Vector{SMat{ZZRingElem}})::SRow{ZZRingElem}
    """
    calculates vector associated with monomial mon
    """
    vec = v0
    degree_mon = degrees(mon)
    for i in length(degree_mon):-1:1
        for j in 1:degree_mon[i]
            vec = mul(vec, transpose(matrices_of_operators[i])) # currently there is no sparse matrix * vector mult
        end
    end
    return vec
end

function highest_calc_sub_monomial(x::Vector{ZZMPolyRingElem}, mon::ZZMPolyRingElem, 
                                    calc_monomials::Dict{ZZMPolyRingElem, Tuple{TVec, Vector{Int}}},
                                    number_of_operators::Int)::ZZMPolyRingElem
    """
    returns the key in calc_monomials that can be extended by the least amount of left-operations to mon
    """
    sub_mon = copy(mon)

    for i in 1:number_of_operators
        while is_divisible_by(sub_mon, x[i])
            if haskey(calc_monomials, sub_mon)
                return sub_mon
            else
                sub_mon /= x[i]
            end
        end
    end
    return sub_mon
end

function calc_new_mon!(x::Vector{ZZMPolyRingElem}, mon::ZZMPolyRingElem, weights::Vector{Vector{Int}}, 
                        matrices_of_operators::Vector{SMat{ZZRingElem}}, number_of_operators::Int,
                        calc_monomials::Dict{ZZMPolyRingElem, Tuple{TVec, Vector{Int}}}, 
                        space::Dict{Vector{Int64}, Oscar.BasisLieHighestWeight.VSBasis}, e::Vector{Vector{Int}}, 
                        cache_size::Int)::SRow{ZZRingElem}
    # calculate vector of mon by extending a previous calculated vector to a
    # monom that differs only by left-multiplication, save results in calc_monomials
    sub_mon = highest_calc_sub_monomial(x, mon, calc_monomials, number_of_operators)
    #println("sub_mon: ", sub_mon)
    sub_mon_cur = copy(sub_mon)
    (vec, weight) = calc_monomials[sub_mon]
    for i in number_of_operators:-1:1
        for k in degrees(sub_mon)[i]:(degrees(mon)[i]-1)
            sub_mon_cur *= x[i]
            weight += weights[i]
            if !haskey(space, weight)
                space[weight] = nullSpace()
            end

            vec = mul(vec, transpose(matrices_of_operators[i])) # currently there is no sparse matrix * vector mult
            if length(calc_monomials) < cache_size
                calc_monomials[sub_mon_cur] = (vec, weight)
            end

            # check if the extended monomial can be deleted from calculated_monomials, i.e. the other possible 
            # extensions by left multiplication with some x[i] are already contained
            can_be_deleted = true
            k = number_of_operators
            for l in 1:number_of_operators
                if degrees(sub_mon_cur - x[i])[l] != 0
                    k = l
                end
            end
            for l in 1:k
                can_be_deleted = can_be_deleted && haskey(calc_monomials, sub_mon_cur - x[i] + x[l])
            end
            if can_be_deleted && sub_mon_cur != x[i]
                delete!(calc_monomials, sub_mon_cur - x[i])
            end
        end
    end
    #println(length(calc_monomials))
    return vec
end

