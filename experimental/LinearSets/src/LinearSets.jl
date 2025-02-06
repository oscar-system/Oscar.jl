module LinearSets

using Oscar


# Need collection of normalized vectors
# needs tweaking if order(base_ring(V)) is not prime....
function _points(V::AbstractAlgebra.FreeModule{T}) where T <: FinFieldElem
    p = order(base_ring(V))
    @req is_prime(p) "Currently only works for V over a prime field"

    B = basis(V)
    n = length(B)

    # TODO : should convert these to linear combinations of B instead of raw conversion.
    [V(digits(p^s + j, base = p, pad = n)) for s in 1:n for j in 0:(p^(s-1)-1)]
end

# want to represent field E as a vector space over a subfield.
function _vector_space(E::FiniteField, phi::Map)
    k = domain(phi)
    @req is_prime(order(k)) "Currently only works for V over a prime field"
    n = degree(E)
    V = vector_space(k,n)
    phi = map_from_func(a -> V(absolute_coordinates(a)), E, k)

    return V, phi
end

# need to store:
#     1. base field F_q so E = F_q^n
#     2. map phi from E -> (Fq)^n
#     3. set S = Points((Fq)^n) as elements of E
# TODO decide: Second argument- subfield k, or an embedding k -> E ?
# NOTE at the moment, we only support F = prime field...
function linear_set_field_attributes!(K::FinField, phi::Map)
    k = domain(phi)
    n = divexact(absolute_degree(K), absolute_degree(k))
    phi = absolute_coordinates
    S =
end function


# nu = first(x for x in E if (x != 0 && !is_square(-x) && is_square(-x-1)))
#
# M = matrix(E, 3, 3 [0 1 0])

end
