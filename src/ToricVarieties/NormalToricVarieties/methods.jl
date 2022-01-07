############################
# Betti numbers
############################

@doc Markdown.doc"""
    betti_number(v::AbstractNormalToricVariety, i::Int)

Compute the `i`-th Betti number of the normal toric variety `v`.
"""
function betti_number(v::AbstractNormalToricVariety, i::Int)
    # check input
    d = dim(v)::Int
    if i > 2*d || i < 0 || isodd(i)
        return fmpz(0)
    end

    # extract vector of currently-known Betti numbers (or create it if necessary)
    betti_numbers = get_attribute!(() -> fill(fmpz(-1),d+1), v, :betti_number)::Vector{fmpz}

    # compute the Betti number if needed
    k = i >> 1 # i is even, so divide by two and use that as index
    if betti_numbers[k+1] == -1
        f_vector::Vector{Int} = pm_object(v).F_VECTOR
        pushfirst!(f_vector, 1)
        betti_numbers[k+1] = fmpz(sum((-1)^(i-k) * binomial(i,k) * f_vector[d - i + 1] for i=k:d))
    end
    
    # return result
    return betti_numbers[k+1]
end
export betti_number
