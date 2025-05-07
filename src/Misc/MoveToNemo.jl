###############################################################################
# A place to accumulate code that should eventually be moved to Nemo.jl
###############################################################################

function minpoly(a::fpFieldElem)
  kx, x = polynomial_ring(parent(a), cached = false)
  return x-a
end

function (F::QQField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{QQPolyRingElem, QQMPolyRingElem}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:QQFieldElem, QQFieldElem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::QQField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{QQFieldElem})
  return F(x.d)
end

function (F::fpField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{fpPolyRingElem, fpMPolyRingElem}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:fpFieldElem, fpFieldElem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::FqField)(x::AbstractAlgebra.Generic.FracFieldElem{T}) where T <: Union{FqPolyRingElem, FqMPolyRingElem}
  num = numerator(x)
  cst_num = constant_coefficient(num)
  denom = denominator(x)
  cst_denom = constant_coefficient(denom)
  if (num != cst_num || denom != cst_denom) throw(InexactError(:FqFieldElem, FqFieldElem, x)) end
  F(cst_num) // F(cst_denom)
end

function (F::fpField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{fpFieldElem})
  return F(x.d)
end

function (F::FqField)(x::AbstractAlgebra.Generic.RationalFunctionFieldElem{FqFieldElem})
  return F(x.d)
end
