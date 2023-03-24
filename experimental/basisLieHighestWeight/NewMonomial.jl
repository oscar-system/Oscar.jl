# main file
#--- bekommt gerade noch ZZ, Short und TVEC aus VectorSpaceBases


# TODO use dimension of weightspace to stop unnecessary calculations
# TODO use UInt8 instead of Int for polynomials again

# module MB2

include("./VectorSpaceBases.jl")
include("./TensorModels.jl")
include("./LieAlgebras.jl")
include("./MonomialOrder.jl")
include("./WeylPolytope.jl")

fromGap = Oscar.GAP.gap_to_julia


function calc_wt(mon, weights)
    """
    
    """
    wt = [0 for i in 1:length(weights[1])]
    for i in 1:length(mon)
        wt .+= mon[i] * weights[i]
    end
    return wt
end

function calc_vec(v0, mon, mats)
    vec = v0
    for i in length(mon):-1:1
        for j in 1:mon[i]
            vec = mul(vec, transpose(mats[i]))
        end
    end
    return vec
end

function highest_calc_sub_monomial(mon::Vector{Int}, calc_monomials)
    """
    returns the key in calc_monomials that can be extended by the least amount of left-operations to mon
    """
    sub_mon = copy(mon)

    m = length(mon)
    for i in 1:m
        while sub_mon[i] > 0
            #println(sub_mon)
            #println(calc_monomials)
            if haskey(calc_monomials, sub_mon)
                return sub_mon
            else
                sub_mon[i] -= 1
            end
        end
    end
    #println(sub_mon)
    return sub_mon # [0 for i in 1:m]
end

function calc_new_mon!(mon, m, wts, mats, calc_monomials, space, e, cache_size)
    # calculate vector of mon by extending a previous calculated vector to a
    # monom that differs only by left-multiplication, save results in calc_monomials
    sub_mon = highest_calc_sub_monomial(mon, calc_monomials)
    #println("sub_mon: ", sub_mon)
    sub_mon_cur = copy(sub_mon)
    if !haskey(calc_monomials, sub_mon_cur)
        println("ERROR IN calc_new_mon")
        println(mon)
        println(sub_mon_cur)
        println(calc_monomials)
    end
    (vec, wt) = calc_monomials[sub_mon]

    # this block cannot be used in MB3 because we don't iterate through monomials in universal order,
    # but instead in the monomial-order for each weightspace
    #if mon == sub_mon # mon already contained so we need to reduce it to 0
    #    return spzeros(ZZ, m), wt
    #end

    #needs_to_be_saved = true
    for i in m:-1:1
        for k in sub_mon[i]:(mon[i]-1)
            sub_mon_cur += e[i]
            wt += wts[i]
            if !haskey(space, wt)
                space[wt] = nullSpace()
            end
            #if isempty(vec.nzind) # v already 0
            #    needs_to_be_saved = false
            #else
            #    vec = mats[i] * vec
            #end
            #vec = addAndReduce!(space[wt], vec)
            #println(sub_mon_cur)
            #if needs_to_be_saved
            #    calc_monomials[sub_mon_cur] = (vec, wt)
            #end
            #mul!(A, vec)
            vec = mats[i] * vec
            if length(calc_monomials) < cache_size
                calc_monomials[sub_mon_cur] = (vec, wt)
            end

            # check if the extended monomial can be deleted from calculated_monomials, i.e. the other possible extensions are already contained
            can_be_deleted = true
            k = m
            for l = 1:m
                if (sub_mon_cur-e[i])[l] != 0
                    k = l
                end
            end
            for l = 1:k
                can_be_deleted = can_be_deleted && haskey(calc_monomials, sub_mon_cur-e[i]+e[l])
            end
            if can_be_deleted && sub_mon_cur != e[i]
                delete!(calc_monomials, sub_mon_cur-e[i])
            end
        end
    end
    # calc_monomials[sub_mon_cur] = (vec, wt) this position is for graded_reverse_lexicographic enough instead of the one above
    return vec, wt
end

