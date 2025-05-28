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
    F = one(U.algebra)
    monomial_set!(F.poly, 1, b)
    for i in 1:length(b)
      mul!(F, inv(q_factorial(b[i], U.qi[i])))
    end

    G = F
    F = sub!(U.bar_automorphism(F), F)
    cf = zero(coefficient_ring(U))
    while !iszero(F)
      exponent_vector!(b, F, length(F))
      elem = _canonical_basis_elem(U, b)
      cf = div!(cf, coeff(F, length(F)), coeff(elem, length(elem)))
      F = submul!(F, elem, cf)

      # cf is a bar-invariant Laurent polynomial, i.e. cf = p(q) / q^k
      # split the coefficient into positive and negative powers of q
      # use postive powers for the canoncial basis element and discard the negative powers
      cf.d.d.num = shift_right!(cf.d.d.num, degree(cf.d.d.den))
      cf.d.d.den = one!(cf.d.d.den)
      G = addmul!(G, elem, cf)
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
