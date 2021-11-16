############################
# Dimensions
############################

@doc Markdown.doc"""
    ith_betti_number(v::AbstractNormalToricVariety, i::Int)

Compute the i-th Betti number of the normal toric variety `v`.
"""
function ith_betti_number(v::AbstractNormalToricVariety, i::Int)
    if isodd(i)
        return 0
    end
    k = div(i, 2)
    f_vector = Vector{Int}(pm_object(v).F_VECTOR)
    pushfirst!(f_vector, 1)
    betti_number = sum((-1)^(i-k) * binomial(i,k) * f_vector[dim(v) - i + 1] for i=k:dim(v))
    return betti_number
    
end
export ith_betti_number
