module AlgClosureFp

using Oscar

import Base: +, -, *, //, hash, show, ==
import Oscar: divexact, add!, mul, mul!, addeq!, sub, data

struct AlgClosure{T} <: AbstractAlgebra.Field
  # T <: FinField
  k::T
  fld::Dict{Int, FinField}
  function AlgClosure(k::T) where T <: FinField
    return new{T}(k, Dict{Int, FinField}(degree(k) => k))
  end
end

function show(io::IO, A::AlgClosure)
  print(io, "Algebraic Closure of $(A.k)")
end

Oscar.base_field(A::AlgClosure) = A.k
Oscar.base_ring(A::AlgClosure) = A.k
Oscar.characteristic(k::AlgClosure) = characteristic(base_field(k))

struct AlgClosureElem{T} <: FieldElem
  # T <: FinFieldElem
  data::FinFieldElem
  parent::AlgClosure{T}
end

Oscar.elem_type(::AlgClosure{T}) where T = AlgClosureElem{T}
Oscar.parent_type(::AlgClosureElem{T}) where T = AlgClosure{T}
Oscar.parent_type(::Type{AlgClosureElem{T}}) where T = AlgClosure{T}

function show(io::IO, a::AlgClosureElem)
  print(io, data(a))
end

function Base.deepcopy_internal(a::AlgClosureElem, d::IdDict)
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


Oscar.zero(A::AlgClosure) = AlgClosureElem(zero(base_field(A)), A)
Oscar.one(A::AlgClosure) = AlgClosureElem(one(base_field(A)), A)

Oscar.parent(a::AlgClosureElem) = a.parent
Oscar.data(a::AlgClosureElem) = a.data

function check_parent(a::AlgClosureElem, b::AlgClosureElem)
  parent(a) == parent(b) || error("incompatible elements")
end

function ext_of_degree(A::AlgClosure, d::Int)
  if haskey(A.fld, d)
    return A.fld[d]
  end
  k = base_ring(A)
  if isa(k, Nemo.GaloisField) || isa(k, FqNmodFiniteField)
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

function Oscar.embed(k::Nemo.GaloisField, K::FqNmodFiniteField)
  @assert characteristic(K) == characteristic(k)
end

+(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(+, a, b), parent(a))
-(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(-, a, b), parent(a))
*(a::AlgClosureElem{T}, b::AlgClosureElem{S}) where S where T = AlgClosureElem(op(*, a, b), parent(a))
//(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(//, a, b), parent(a))
divexact(a::AlgClosureElem, b::AlgClosureElem) = AlgClosureElem(op(divexact, a, b), parent(a))
==(a::AlgClosureElem, b::AlgClosureElem) = op(==, a, b)
Oscar.iszero(a::AlgClosureElem) = iszero(data(a))
Oscar.isone(a::AlgClosureElem) = isone(data(a))
-(a::AlgClosureElem) = AlgClosureElem(-data(a), parent(a))

function Oscar.roots(a::AlgClosureElem, b::Int)
  ad = data(a)
  kx, x = PolynomialRing(parent(ad), cached = false)
  f = x^b-ad
  lf = factor(f)
  d = mapreduce(degree, lcm, keys(lf.fac), init = 1)
  d = lcm(d, degree(parent(ad)))
  K = ext_of_degree(parent(a), d)
  r = roots(f, K)
  return [AlgClosureElem(x, parent(a)) for x = r]
end

function Oscar.roots(a::Generic.Poly{AlgClosureElem{T}}) where T
  A = base_ring(a)
  b = minimize(FinField, collect(coefficients(a)))
  kx, x = PolynomialRing(parent(b[1]), cached = false)
  f = kx(b)
  lf = factor(f)
  d = mapreduce(degree, lcm, keys(lf.fac), init = 1)
  d = lcm(d, degree(parent(b[1])))
  K = ext_of_degree(A, d)
  r = roots(f, K)
  return [AlgClosureElem(x, A) for x = r]
end


function Oscar.minpoly(a::AlgClosureElem)
  return minpoly(data(a))
end

function Oscar.minpoly(a::gfp_elem)
  kx, x = PolynomialRing(parent(a), cached = false)
  return x-a
end

function Oscar.degree(a::AlgClosureElem)
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

function Base.hash(a::AlgClosureElem, u::UInt)
  b = minimize(a)
  return hash(data(b), u)
end

function Oscar.gmodule(::Type{FinField}, C::GModule{<:Any, <:Generic.FreeModule{<:AlgClosureElem{<:FinField}}})

  d = dim(C)
  l = 1
  for g = C.ac
    l = lcm(l, lcm(collect(map_entries(x->Hecke.degree(parent(x.data)), mat(g)))))
  end
  K = ext_of_degree(base_ring(C), l)
  return gmodule(K, C)
end

function Oscar.map_entries(K::FinField, M::MatElem{<:AlgClosureElem})
  N = zero_matrix(K, nrows(M), ncols(M))
  for i=1:nrows(M)
    for j=1:ncols(M)
      embed(parent(data(M[i,j])), K)
      N[i,j] = K(data(M[i,j]))
    end
  end
  return N
end

function Oscar.gmodule(K::FinField, C::GModule{<:Any, <:Generic.FreeModule{<:AlgClosureElem{<:FinField}}})

  d = dim(C)
  F = free_module(K, d)
  if d == 0 
    h = hom(F, F, elem_type(F)[])
    return gmodule(F, group(C), typeof(h)[hom(F, F, map_entries(K, mat(x))) for x = C.ac])
  end
  return gmodule(F, group(C), [hom(F, F, map_entries(K, mat(x))) for x = C.ac])
end

function Oscar.gmodule(K::FinField, C::GModule{<:Any, <:Generic.FreeModule{<:FinFieldElem}})

  d = dim(C)
  F = free_module(K, d)
  if d == 0 
    h = hom(F, F, elem_type(F)[])
    return gmodule(F, group(C), typeof(h)[hom(F, F, map_entries(K, mat(x))) for x = C.ac])
  end
  return gmodule(F, group(C), [hom(F, F, map_entries(K, mat(x))) for x = C.ac])
end


function Oscar.GModuleFromGap.hom_base(C::T, D::T) where T <: GModule{<:Any, <:Generic.FreeModule{<:AlgClosureElem{<:FinField}}}

  C1 = gmodule(FinField, C)
  D1 = gmodule(FinField, D)
  Cf = degree(base_ring(C1))
  Df = degree(base_ring(D1))
  l = lcm(Cf, Df)
  K = ext_of_degree(base_ring(C), l)
  if l != Cf
    C1 = gmodule(K, C1)
  end
  if l != Df
    D1 = gmodule(K, D1)
  end
  h = Oscar.GModuleFromGap.hom_base(C1, D1)
  if length(h) == 0
    return h
  end
  return map(x->map_entries(base_ring(C), x), h)
end


end # AlgClosureFp
