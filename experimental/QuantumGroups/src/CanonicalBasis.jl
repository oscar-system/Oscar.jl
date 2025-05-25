@doc raw"""
    canonical_basis_elem(U::QuantumGroup, b::Vector{Int}) -> QuantumGroupElem
    
Return the element of the canonical basis in the negative part of `U`
corresponding to the Lusztig datum `b`.

# Examples
```jldoctest
julia> U = quantum_group(:A, 2);

julia> canonical_basis_elem(U, [0, 1, 0])
q^-1*F[1]*F[3] + F[2]
"""
function canonical_basis_elem(U::QuantumGroup, b::Vector{Int})
  return QuantumGroupElem(U, _canonical_basis_elem(U, b.ref.mem))
end

function _canonical_basis_elem(U::QuantumGroup, b::Memory{Int})
  return get!(U.canonical_basis, b) do
    bar = bar_involution(U)
    F = zero(U.algebra)
    add_monomial!(F.poly, b)
    for i in 1:length(b)
      mul!(F, inv(quantum_factorial(b[i], U.qi[i])))
    end

    G = F
    F = sub!(bar(F), F)
    while !iszero(F)
      exponent_vector!(b, F, length(F))
      elem = _canonical_basis_elem(U, b)
      cf1 = numerator(div(coeff(F, length(F)), coeff(elem, length(elem))))

      # split the coefficient into positive and negative powers of q
      # use postive powers for the canoncial basis element and discard the negative powers
      # this corresponds to b and bar(b)
      cf2 = coefficient_ring(U)(
        parent(cf1.poly)([coeff(cf1, n) for n in 0:degree(cf1.poly)])
      )

      G = addmul!(G, elem, cf2)
      F = submul!(F, elem, coefficient_ring(U)(cf1))
    end

    return G
  end
end

@doc raw"""
    canonical_basis_expansion(x::QuantumGroupElem) -> Vector{Tuple{QuantumFieldElem,Vector{Int}}}
    
Return the expansion into the canonical basis of an element `x` lying in the negative part `parent(x)`.
The expansion is given as a vector of tuples,
where each tuple contains the coefficient and the Lusztig datum of a canonical basis element.
"""
function canonical_basis_expansion(x::QuantumGroupElem)
  U = parent(x)
  rep = Tuple{QuantumFieldElem,Vector{Int}}[]

  y = deepcopy(x.elem)
  b = Memory{Int}(undef, ngens(U))
  while !iszero(y)
    exponent_vector!(b, y, length(y))
    elem = _canonical_basis_elem(U, b)

    coeff = coeff(y, length(y)) / coeff(elem, length(elem))
    push!(rep, (coeff, exp))
    y = submul!(y, elem, coeff)
  end

  return rep
end
