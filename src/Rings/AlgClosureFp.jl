###############################################################################
# 
#  Algebraic closure of finite fields
#
###############################################################################

# This is an implementation of the algebraic closure of finite fields,
# which is modeled as the union of finite fields.

module AlgClosureFp

using ..Oscar

import Base: +, -, *, //, ==, deepcopy_internal, hash, isone, iszero, one,
  parent, show, zero

import ..Oscar: pretty, Lowercase

import ..Oscar: algebraic_closure, base_field, base_ring, base_ring_type, characteristic, data, degree, divexact,
  elem_type, embedding, has_preimage_with_preimage, IntegerUnion, is_unit, map_entries,
  minpoly, parent_type, promote_rule, roots

"""
    AlgClosure{T} <: AbstractAlgebra.Field

Type for the algebraic closure of a finite field.

See [`algebraic_closure`](@ref).
"""
struct AlgClosure{T} <: AbstractAlgebra.Field
  # T <: FinField
  k::T
  fld::Dict{Int, FinField} # Cache for the finite fields
  function AlgClosure(k::T) where T <: FinField
    return new{T}(k, Dict{Int, FinField}(degree(k) => k))
  end
end

function show(io::IO, A::AlgClosure)
  io = pretty(io)
  print(io, "Algebraic closure of ", Lowercase(), A.k)
end

base_field(A::AlgClosure) = A.k
base_ring(A::AlgClosure) = A.k
base_ring_type(::Type{AlgClosure{T}}) where {T} = T
characteristic(k::AlgClosure) = characteristic(base_field(k))

struct AlgClosureElem{T} <: FieldElem
  # T <: FinField
  data::FinFieldElem
  parent::AlgClosure{T}
end

elem_type(::Type{AlgClosure{T}}) where T = AlgClosureElem{T}
parent_type(::Type{AlgClosureElem{T}}) where T = AlgClosure{T}

Oscar.canonical_unit(a::AlgClosureElem) = is_zero(a) ? one(a) : a

function show(io::IO, a::AlgClosureElem)
  print(io, data(a))
end

function deepcopy_internal(a::AlgClosureElem, d::IdDict)
  return AlgClosureElem(deepcopy_internal(data(a), d), parent(a))
end

(A::AlgClosure)(a::IntegerUnion) = AlgClosureElem(A.k(a), A)
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

#TODO: Guarantee to return a field of the same type as `base_ring(A)`?
# (Then `fpField` cannot be supported as `base_ring(A)`)
@doc raw"""
    ext_of_degree(A::AlgClosure, d::Int)

Return a finite field `F` of order `p^d`
where `p` is the characteristic of `K`.
This field is compatible with `A` in the sense that `A(x)` returns the
element of `A` that corresponds to the element `x` of `F`.
"""
function ext_of_degree(A::AlgClosure, d::Int)
  if haskey(A.fld, d)
    return A.fld[d]
  end

  #see if the subfield lattice already has the field...
  for x = values(A.fld)
    degree(x) % d == 0 || continue
    for (deg, mK) = Nemo.subfields(x)
      if d == deg
        A.fld[d] = K = domain(mK[1])
        return K
      end
    end
  end
    
  k = base_ring(A)
  if isa(k, fpField) || isa(k, fqPolyRepField)
    K = Nemo.Native.GF(Int(characteristic(k)), d, cached = false)
  elseif isa(k, FqField)
    K = GF(characteristic(k), d, cached = false)
  else
    K = Nemo.Native.GF(characteristic(k), d, cached = false)
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

+(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(+, a, b), parent(a))
-(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(-, a, b), parent(a))
#TODO: do we really want to support different types here? (implies different parents)
*(a::AlgClosureElem{T}, b::AlgClosureElem{S}) where S where T = AlgClosureElem(op(*, a, b), parent(a))
//(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(//, a, b), parent(a))
divexact(a::AlgClosureElem, b::AlgClosureElem; check = true) = AlgClosureElem(op(divexact, a, b), parent(a))
==(a::AlgClosureElem, b::AlgClosureElem) = op(==, a, b)
iszero(a::AlgClosureElem) = iszero(data(a))
isone(a::AlgClosureElem) = isone(data(a))
-(a::AlgClosureElem) = AlgClosureElem(-data(a), parent(a))


################################################################################
#
#  Ad hoc binary operations
#
################################################################################

*(a::IntegerUnion, b::AlgClosureElem) = AlgClosureElem(data(b)*a, parent(b))

*(a::AlgClosureElem, b::IntegerUnion) = b*a

+(a::IntegerUnion, b::AlgClosureElem) = AlgClosureElem(data(b) + a, parent(b))

+(a::AlgClosureElem, b::IntegerUnion) = b + a

-(a::IntegerUnion, b::AlgClosureElem) = AlgClosureElem(-(a, data(b)), parent(b))

-(a::AlgClosureElem, b::IntegerUnion) = AlgClosureElem(-(data(a), b), parent(a))

//(a::AlgClosureElem, b::Integer) = AlgClosureElem(data(a)//b, parent(a))
//(a::AlgClosureElem, b::ZZRingElem) = AlgClosureElem(data(a)//b, parent(a))


is_unit(a::AlgClosureElem) = !iszero(a)


###############################################################################
#
#   Functions for computing roots
#
###############################################################################

function roots(a::AlgClosureElem, b::Int)
  ad = data(a)
  kx, x = polynomial_ring(parent(ad); cached = false)
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
  kx, x = polynomial_ring(parent(b[1]); cached = false)
  f = kx(b)
  lf = factor(f)
  d = mapreduce(degree, lcm, keys(lf.fac), init = 1)
  d = lcm(d, degree(parent(b[1])))
  K = ext_of_degree(A, d)
  r = roots(K, f)
  return [AlgClosureElem(x, A) for x = r]
end


#TODO: Does this really make sense?
# K = algebraic_closure(GF(3, 1))
# F2 = ext_of_degree(K, 2)
# a = gen(F2); fa = minpoly(a); fa(a)  # works
# c = K(a); fc = minpoly(c); fc(c)  # does not work
function minpoly(a::AlgClosureElem)
  return minpoly(data(a))
end

# Note: We want the degree of the smallest finite field that contains `a`.
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
  @assert allequal(parent, a)
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


################################################################################
#
#  Creation of the field
#
################################################################################

@doc raw"""
    algebraic_closure(F::FinField)

Let `F` be a prime field of order `p`.
Return a field `K` that is the union of finite fields of order `p^d`,
for all positive integers `d`.
The degree `d` extension of `F` can be obtained as `ext_of_degree(K, d)`.

`K` is cached in `F`, and the fields returned by `ext_of_degree` are
cached in `K`.

# Examples
```jldoctest; setup = :(using Oscar)
julia> K = algebraic_closure(GF(3));

julia> F2 = ext_of_degree(K, 2);

julia> F3 = ext_of_degree(K, 3);

julia> x = K(gen(F2)) + K(gen(F3));

julia> degree(x)
6
```
"""
@attr AlgClosure{T} function algebraic_closure(F::T) where T <: FinField
  @req is_prime(order(F)) "only for finite prime fields"
  return AlgClosure(F)
end

function embedding(k::T, K::AlgClosure{T}) where T <: FinField
  @req characteristic(k) == characteristic(K) "incompatible characteristics"
  return MapFromFunc(k, K, K, k)
end

function has_preimage_with_preimage(mp::MapFromFunc{T, AlgClosure{S}}, elm::AlgClosureElem{S}) where T <: FinField where S <: FinField
  F = domain(mp)
  mod(degree(F), degree(elm)) != 0 && return false, zero(F)
  return true, preimage(mp, elm)
end

### Conformance test element generation
function ConformanceTests.generate_element(K::AlgClosure{T}) where T <: FinField
  d = rand(1:8)
  F = ext_of_degree(K, d)
  return K(rand(F))
end

end # AlgClosureFp

import .AlgClosureFp:
       AlgClosure,
       AlgClosureElem,
       ext_of_degree

export AlgClosure,
       AlgClosureElem,
       ext_of_degree
