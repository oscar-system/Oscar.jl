############################
# Betti numbers
############################

@doc Markdown.doc"""
    betti_number(v::AbstractNormalToricVariety, i::Int)

Compute the i-th Betti number of the normal toric variety `v`.
"""
function betti_number(v::AbstractNormalToricVariety, i::Int)
    # check input
    if i > 2*dim(v) || i < 0
        return 0
    end
    
    # extract vector of currently-known Betti numbers (or create it if necessary)
    if !has_attribute(v, :betti_number)
        betti_numbers = fill(fmpz(-1),2*dim(v)+1)
    else
        betti_numbers = get_attribute(v, :betti_number)::Vector{fmpz}
    end
    
    # compute the Betti number if needed
    if betti_numbers[i+1] == -1
        if isodd(i)
            betti_numbers[i+1] = fmpz(0)
        else
            k = div(i, 2)
            f_vector = Vector{Int}(pm_object(v).F_VECTOR)
            pushfirst!(f_vector, 1)
            betti_numbers[i+1] = fmpz(sum((-1)^(i-k) * binomial(i,k) * f_vector[dim(v) - i + 1] for i=k:dim(v)))
        end
        set_attribute!(v, :betti_number, betti_numbers)
    end
    
    # return result
    return get_attribute(v, :betti_number)[i+1]
end
export betti_number
