function tor(M1::SparseFPModule{T}, M2::SparseFPModule{T}) where {U<:MPolyComplementOfPrimeIdeal, T<:MPolyLocRingElem{<:Any, <:Any, <:Any, <:Any, U}}
  f1, _ = free_resolution(SimpleFreeResolution, M1)
  M1oM2 = tensor_product(free_resolution(SimpleFreeResolution, M1)[1], ZeroDimensionalComplex(M2))
end

function homology(hc::AbsHyperComplex, i::Int)
  dim(hc) == 1 || error("dimension of the hypercomplex must be one")
  return homology(hc, 1, (i,))
end

function kernel(hc::AbsHyperComplex, i::Int)
  dim(hc) == 1 || error("dimension of the hypercomplex must be one")
  return kernel(hc, 1, (i,))
end

function boundary(hc::AbsHyperComplex, i::Int)
  dim(hc) == 1 || error("dimension of the hypercomplex must be one")
  return boundary(hc, 1, (i,))
end

function homology(hc::AbsHyperComplex, p::Int, i::Tuple)
  return homology(underlying_complex(hc), p, i)
end

function kernel(hc::AbsHyperComplex, p::Int, i::Tuple)
  return kernel(underlying_complex(hc), p, i)
end

function boundary(hc::AbsHyperComplex, p::Int, i::Tuple)
  return boundary(underlying_complex(hc), p, i)
end

function tor(M1::SparseFPModule{T}, M2::SparseFPModule{T}, i::Int) where {U<:MPolyComplementOfPrimeIdeal, T<:MPolyLocRingElem{<:Any, <:Any, <:Any, <:Any, U}}
  return homology(tor(M1, M2), i)
end

function euler_characteristic(hc::AbsHyperComplex)
  dim(hc) == 1 || error("complex must be one-dimensional")
  has_lower_bound(hc) && has_upper_bound(hc) || error("complex must be bounded")
  return sum((-1)^i*length(homology(hc, i)) for i in lower_bound(hc):upper_bound(hc); init=0)
end

