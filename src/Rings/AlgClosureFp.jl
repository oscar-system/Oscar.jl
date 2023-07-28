###############################################################################
# 
#  Algebraic closure of finite fields
#
###############################################################################

# This is an implementation of the algebraic closure of finite fields,
# which is modelled as the union of finite fields.

module AlgClosureFp

using ..Oscar

import Base: +, -, *, //, ==, deepcopy_internal, hash, isone, iszero, one,
  parent, show, zero

import ..Oscar: base_field, base_ring, characteristic, data, degree, divexact,
  elem_type, map_entries, minpoly, parent_type, roots

struct AlgClosure{T} <: AbstractAlgebra.Field
  # T <: FinField
  k::T
  fld::Dict{Int, FinField} # Cache for the finite fields
  function AlgClosure(k::T) where T <: FinField
    return new{T}(k, Dict{Int, FinField}(degree(k) => k))
  end
end

function show(io::IO, A::AlgClosure)
  print(io, "Algebraic Closure of $(A.k)")
end

base_field(A::AlgClosure) = A.k
base_ring(A::AlgClosure) = A.k
characteristic(k::AlgClosure) = characteristic(base_field(k))

struct AlgClosureElem{T} <: FieldElem
  # T <: FinFieldElem
  data::FinFieldElem
  parent::AlgClosure{T}
end

elem_type(::Type{AlgClosure{T}}) where T = AlgClosureElem{T}
elem_type(::AlgClosure{T}) where T = AlgClosureElem{T}
parent_type(::AlgClosureElem{T}) where T = AlgClosure{T}
parent_type(::Type{AlgClosureElem{T}}) where T = AlgClosure{T}

function show(io::IO, a::AlgClosureElem)
  print(io, data(a))
end

function deepcopy_internal(a::AlgClosureElem, d::IdDict)
  return AlgClosureElem(data(a), parent(a))
end

(A::AlgClosure)(a::Int) = AlgClosureElem(A.k(a), A)
(A::AlgClosure)(a::AlgClosureElem) = a
(A::AlgClosure)() = A(0)
function (A::AlgClosure)(a::FinFieldElem)
  @assert characteristic(parent(a)) == characteristic(A)
  if haskey(A.fld, degree(parent(a)))
    @assert A.fld[degree(parent(a))] == parent(a)
  end
  return AlgClosureElem(a, A)
end


zero(A::AlgClosure) = AlgClosureElem(zero(base_field(A)), A)
one(A::AlgClosure) = AlgClosureElem(one(base_field(A)), A)

parent(a::AlgClosureElem) = a.parent
data(a::AlgClosureElem) = a.data

function check_parent(a::AlgClosureElem, b::AlgClosureElem)
  parent(a) == parent(b) || error("incompatible elements")
end

function ext_of_degree(A::AlgClosure, d::Int)
  if haskey(A.fld, d)
    return A.fld[d]
  end
  k = base_ring(A)
  if isa(k, Nemo.fpField) || isa(k, fqPolyRepField)
    K = GF(Int(characteristic(k)), d, cached = false)
  else
    K = GF(characteristic(k), d, cached = false)
  end
  A.fld[d] = K
  return K
end

function op(f::Function, a::AlgClosureElem, b::AlgClosureElem)
  check_parent(a, b)
  ad = data(a)
  bd = data(b)
  if parent(ad) == parent(bd)
    return f(ad,bd)
  end

  l = lcm(degree(parent(ad)), degree(parent(bd)))
  k = ext_of_degree(parent(a), l)
  embed(parent(ad), k)
  embed(parent(bd), k)
  return f(k(ad), k(bd))
end

#T the following belongs to Nemo and should be moved there
function Oscar.embed(k::Nemo.fpField, K::fqPolyRepField)
  @assert characteristic(K) == characteristic(k)
end

+(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(+, a, b), parent(a))
-(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(-, a, b), parent(a))
*(a::AlgClosureElem{T}, b::AlgClosureElem{S}) where S where T = AlgClosureElem(op(*, a, b), parent(a))
//(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(//, a, b), parent(a))
divexact(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(divexact, a, b), parent(a))
==(a::AlgClosureElem, b::AlgClosureElem) = op(==, a, b)
iszero(a::AlgClosureElem) = iszero(data(a))
isone(a::AlgClosureElem) = isone(data(a))
-(a::AlgClosureElem) = AlgClosureElem(-data(a), parent(a))

function roots(a::AlgClosureElem, b::Int)
  ad = data(a)
  kx, x = polynomial_ring(parent(ad), cached = false)
  f = x^b-ad
  lf = factor(f)
  d = mapreduce(degree, lcm, keys(lf.fac), init = 1)
  d = lcm(d, degree(parent(ad)))
  K = ext_of_degree(parent(a), d)
  r = roots(K, f)
  return [AlgClosureElem(x, parent(a)) for x = r]
end

function roots(a::Generic.Poly{AlgClosureElem{T}}) where T
  A = base_ring(a)
  b = minimize(FinField, collect(coefficients(a)))
  kx, x = polynomial_ring(parent(b[1]), cached = false)
  f = kx(b)
  lf = factor(f)
  d = mapreduce(degree, lcm, keys(lf.fac), init = 1)
  d = lcm(d, degree(parent(b[1])))
  K = ext_of_degree(A, d)
  r = roots(K, f)
  return [AlgClosureElem(x, A) for x = r]
end


function minpoly(a::AlgClosureElem)
  return minpoly(data(a))
end

function minpoly(a::fpFieldElem)
  kx, x = polynomial_ring(parent(a), cached = false)
  return x-a
end

function degree(a::AlgClosureElem)
  #TODO: via Frobenius? as a fixed s.th.?
  return degree(minpoly(data(a)))
end

function minimize(a::AlgClosureElem)
  f = minpoly(a)
  k = ext_of_degree(parent(a), degree(f))
  embed(k, parent(data(a)))
  return AlgClosureElem(k(data(a)), parent(a))
end

function minimize(::Type{FinField}, a::AlgClosureElem)
  return data(minimize(a))
end

function minimize(::Type{FinField}, a::AbstractArray{<:AlgClosureElem})
  if length(a) == 0
    return a
  end
  @assert all(x->parent(x) == parent(a[1]), a)
  da = map(degree, a)
  l = reduce(lcm, da)
  k = ext_of_degree(parent(a[1]), l)
  b = elem_type(k)[]
  for i = eachindex(a)
    if da[i] < l
      embed(parent(data(a[i])), k)
      push!(b, k(data(a[i])))
    elseif da[i] == l
      push!(b, data(a[i]))
    else
      embed(k, parent(data(a[i])))
      push!(b, k(data(a[i])))
    end
  end
  return b
end

function (F::FinField)(a::AlgClosureElem)
  b = minimize(FinField, a)
  embed(parent(b), F)
  return F(b)
end

function hash(a::AlgClosureElem, u::UInt)
  b = minimize(a)
  return hash(data(b), u)
end

function map_entries(K::FinField, M::MatElem{<:AlgClosureElem})
  N = zero_matrix(K, nrows(M), ncols(M))
  for i=1:nrows(M)
    for j=1:ncols(M)
      embed(parent(data(M[i,j])), K)
      N[i,j] = K(data(M[i,j]))
    end
  end
  return N
end

end # AlgClosureFp

import .AlgClosureFp:
       AlgClosure,
       AlgClosureElem

export AlgClosure,
       AlgClosureElem
