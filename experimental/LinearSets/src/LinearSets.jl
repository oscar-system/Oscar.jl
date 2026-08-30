module LinearSets

using Oscar

# should use something other than defining_polynomial when v is not over the prime field..
function (k::FqField)(v::AbstractAlgebra.Generic.FreeModuleElem)
  return k(parent(defining_polynomial(k))( Array(v.v)[1:degree(k)] ))
end

# Gives collection of normalized vectors corresponding to 1-dimensional subspaces
# TODO needs tweaking if order(base_ring(V)) is not prime....
function _points(V::AbstractAlgebra.Generic.FreeModule{<:FinFieldElem})
    k = base_ring(V)
    p = order(k)
    @req is_prime(p) "Currently only works for V over a prime field"

    B = basis(V)
    n = length(B)

    # TODO : should convert these to linear combinations of B instead of raw conversion.
    # But this creates overhead, maybe only do this if B is nonstandard...
    # TODO : is it faster to concatenate with the prefix [0...1], rather then compute p^s?
    [V(reverse(digits(p^s + j, base = Int(p), pad = n))) for s in 0:(n-1) for j in 0:(p^s - 1)]
end

# want to represent field E as a vector space over a subfield.
function _vector_space(E::FinField, phi::Map)
    k = domain(phi)
    @req is_prime(order(k)) "Currently only works for V over a prime field"
    n = degree(E)
    V = vector_space(k,n)
    rho = (a -> V(absolute_coordinates(a)))
    rhoI = (v -> E(v))

    return V, MapFromFunc(E, V, rho, rhoI)
end

# need to store:
#     1. base field F_q so E = F_q^n
#     2. map phi from E -> (Fq)^n
#     3. set S = Points((Fq)^n) as elements of E
# TODO decide: Second argument- subfield k, or an embedding k -> E ?
# NOTE at the moment, we only support F = prime field...
# function linear_set_field_attributes!(E::FinField, phi::Map)
#     k = domain(phi)
#     V, rho = _vector_space(E, phi)
#     n = divexact(absolute_degree(K), absolute_degree(k))
#     S = _points(V)
# end


# nu = first(x for x in E if (x != 0 && !is_square(-x) && is_square(-x-1)))
#
# M = matrix(E, 3, 3 [0 1 0])

end
