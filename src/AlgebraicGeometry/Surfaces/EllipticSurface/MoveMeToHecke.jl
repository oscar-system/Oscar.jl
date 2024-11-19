
function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:MPolyRingElem}, a::Vector{T}) where {T<:RingElem}
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end

function evaluate(f::AbstractAlgebra.Generic.FracFieldElem{<:PolyRingElem}, a::RingElem)
  return evaluate(numerator(f), a)//evaluate(denominator(f), a)
end

function extend_domain_to_fraction_field(phi::Map{<:MPolyRing, <:Ring})
  ext_dom = fraction_field(domain(phi))
  return MapFromFunc(ext_dom, codomain(phi), x->phi(numerator(x))*inv(phi(denominator(x))))
end


function iszero(P::EllipticCurvePoint)
  return iszero(P[1]) && isone(P[2]) && iszero(P[3])
end

@doc raw"""
    extended_ade(ADE::Symbol, n::Int)

Return the dual intersection matrix of an extended ade Dynkin diagram
as well as the isotropic vector (with positive coefficients in the roots).
"""
function extended_ade(ADE::Symbol, n::Int)
  R = change_base_ring(ZZ,gram_matrix(root_lattice(ADE,n)))
  G = block_diagonal_matrix([ZZ[2;],R])
  if ADE == :E && n == 8
    G[1,n] = -1
    G[n,1] = -1
  end
  if ADE == :E && n == 7
    G[1,2] = -1
    G[2,1] = -1
  end
  if ADE == :E && n == 6
    G[1,n+1] = -1
    G[n+1,1] = -1
  end
  if ADE == :A && n > 0
    G[1,2] = -1
    G[2,1] = -1
    G[1,n+1] = -1
    G[n+1,1] = -1
  end
  if ADE == :A && n ==1 0
    G[1,2]= -2
    G[2,1] = -2
  end
  if ADE == :D
    G[1,n] = -1
    G[n,1] = -1
  end
  @assert rank(G) == n
  return -G, kernel(G; side = :left)
end
