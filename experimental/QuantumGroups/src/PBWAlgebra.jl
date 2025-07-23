###############################################################################
#
#   Constructors
#
###############################################################################

function pbw_algebra(
  R::AbstractAlgebra.MPolyRing{T}, rels::Vector{<:AbstractAlgebra.MPolyRingElem{T}}
) where {T<:FieldElem}
  S = MPolyRing{T}(ngens(R), coefficient_ring(R), false, symbols(R))
  rels2 = Vector{MPolyRingElem{T}}(undef, length(rels))
  for i in 1:length(rels)
    rels2[i] = zero(S)
    for j in 1:length(rels[i])
      rels2[i] = add!(rels2[i], S(coeff(rels[i], j), exponent_vector(rels[i], j)))
    end
  end

  return PBWAlgebra(S, rels2)
end

function pbw_algebra(
  R::AbstractAlgebra.MPolyRing{T}, rels::Matrix{<:AbstractAlgebra.MPolyRingElem{T}}
) where {T<:FieldElem}
  N = ngens(R)
  return pbw_algebra(R, [rels[i, j] for i in 1:N for j in (i + 1):N])
end

###############################################################################
#
#   Type system
#
###############################################################################

function elem_type(::Type{PBWAlgebra{T}}) where {T}
  return PBWAlgebraElem{T}
end

function parent_type(::Type{PBWAlgebraElem{T}}) where {T}
  return PBWAlgebra{T}
end

function is_domain_type(::Type{PBWAlgebraElem{T}}) where {T}
  return false
end

function is_exact_type(::Type{PBWAlgebraElem{T}}) where {T}
  return is_exact_type(T)
end

###############################################################################
#
#   Accessors
#
###############################################################################

function coefficient_ring(A::PBWAlgebra)
  return coefficient_ring(A.R)
end

function parent(x::PBWAlgebraElem)
  return x.parent
end

###############################################################################
#
#   Basic
#
###############################################################################

function (A::PBWAlgebra)()
  return zero(A)
end

function (A::PBWAlgebra)(a::Integer)
  return PBWAlgebraElem(A, A.R(a))
end

function (A::PBWAlgebra{T})(a::T) where {T}
  return PBWAlgebraElem(A, A.R(a))
end

function (A::PBWAlgebra{T})(x::PBWAlgebraElem{T}) where {T}
  @req A === parent(x) "parent mismatch"
  return x
end

function gen(A::PBWAlgebra, i::Int)
  return PBWAlgebraElem(A, gen(A.R, i))
end

function gens(A::PBWAlgebra)
  return [gen(A, i) for i in 1:ngens(A)]
end

function ngens(A::PBWAlgebra)
  return ngens(A.R)
end

function set!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  check_parent(x, y)
  x.poly = set!(x.poly, y.poly)
  return x
end

function swap!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  x.poly, y.poly = y.poly, x.poly
  return x
end

###############################################################################
#
#   Copying / Hashing
#
###############################################################################

function deepcopy_internal(A::PBWAlgebra, dict::IdDict)
  return get!(dict, A) do
    PBWAlgebra(deepcopy_internal(A.R, dict), deepcopy_internal(A.mult, dict))
  end
end

function deepcopy_internal(x::PBWAlgebraElem{T}, dict::IdDict) where {T}
  return get!(dict, x) do
    PBWAlgebraElem(parent(x), deepcopy_internal(x.poly, dict)::MPolyRingElem{T})
  end
end

function hash(A::PBWAlgebra, h::UInt)
  b = 0x1cd6ee4308df324a % UInt
  h = hash(A.R, h)
  # h = hash(A.mult, h) hashing does not work, due to undef entries

  return xor(h, b)
end

function hash(x::PBWAlgebraElem, h::UInt)
  b = 0xa080ada44dcea378 % UInt
  h = hash(x.parent, h)
  h = hash(x.poly, h)

  return xor(h, b)
end

###############################################################################
#
#   Printing
#
###############################################################################

function Base.show(io::IO, A::PBWAlgebra)
  print(io, "PBW algebra over ", coefficient_ring(A))
end

function Base.show(io::IO, x::PBWAlgebraElem)
  print(io, x.poly)
end

###############################################################################
#
#   Coefficients
#
###############################################################################

function coeff(x::PBWAlgebraElem, i::Int)
  return coeff(x.poly, i)
end

function setcoeff!(x::PBWAlgebraElem{T}, i::Int, a::T) where {T}
  return setcoeff!(x.poly, i, a)
end

###############################################################################
#
#   Arithemtic
#
###############################################################################

function one(A::PBWAlgebra)
  return PBWAlgebraElem(A, one(A.R))
end

function one(x::PBWAlgebraElem)
  return one(parent(x))
end

function one!(x::PBWAlgebraElem)
  x.poly = one!(x.poly)
  return x
end

function zero(A::PBWAlgebra)
  return PBWAlgebraElem(A, zero(A.R))
end

function zero(x::PBWAlgebraElem)
  return zero(parent(x))
end

function zero!(x::PBWAlgebraElem)
  x.poly = zero!(x.poly)
  return x
end

function length(x::PBWAlgebraElem)
  return length(x.poly)
end

function exponent(x::PBWAlgebraElem, i::Int, j::Int)
  return exponent(x.poly, i, j)
end

function set_exponent!(x::PBWAlgebraElem, i::Int, j::Int, e::Int)
  return set_exponent!(x.poly, i, j, e)
end

function exponent_vector!(exp::Vector{Int}, x::PBWAlgebraElem, i::Int)
  return exponent_vector!(exp, x.poly, i)
end

function exponent_vector(x::PBWAlgebraElem, i::Int)
  return exponent_vector(x.poly, i)
end

function set_exponent_vector!(x::PBWAlgebraElem, i::Int, exp::Vector{Int})
  set_exponent_vector!(x.poly, i, exp)
  return x
end

function leading_monomial(x::PBWAlgebraElem)
  return PBWAlgebraElem(parent(x), leading_monomial(x.poly))
end

###############################################################################
#
#   Comparison
#
###############################################################################

function Base.:(==)(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  check_parent(x, y)
  return x.poly == y.poly
end

function iszero(x::PBWAlgebraElem)
  return iszero(x.poly)
end

function isone(x::PBWAlgebraElem)
  return isone(x.poly)
end

function is_unit(x::PBWAlgebraElem)
  return is_unit(x.poly)
end

function is_gen(x::PBWAlgebraElem)
  return is_gen(x.poly)
end

function is_gen_with_index(x::PBWAlgebraElem)
  return is_gen_with_index(x.poly)
end

###############################################################################
#
#   Arithemtic
#
###############################################################################

function Base.:*(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  check_parent(x, y)
  return mul!(zero(x), x, y)
end

function Base.:*(x::PBWAlgebraElem{T}, a::T) where {T<:FieldElem}
  return mul!(zero(x), x, a)
end

function Base.:*(a::T, x::PBWAlgebraElem{T}) where {T<:FieldElem}
  return mul!(zero(x), x, a)
end

function Base.:+(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  check_parent(x, y)
  return add!(zero(x), x, y)
end

function Base.:-(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T<:FieldElem}
  check_parent(x, y)
  return sub!(zero(x), x, y)
end

function Base.:-(x::PBWAlgebraElem)
  return neg!(zero(x), x)
end

function Base.:^(x::PBWAlgebraElem, n::Int)
  return pow!(zero(x), x, n)
end

function divexact_right(x::PBWAlgebraElem{T}, a::T) where {T}
  return PBWAlgebraElem(parent(x), divexact(x.poly, a))
end

function divexact_left(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  error("not implemented")
end

function divexact_right(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  error("not implemented")
end

###############################################################################
#
#   Unsafe arithmetic
#
###############################################################################

function add!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  z.poly = add!(z.poly, x.poly, y.poly)
  return z
end

function addmul!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}, a::T) where {T}
  x.poly = addmul!(x.poly, y.poly, a)
  return x
end

function div!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, a::T) where {T}
  z.poly = div!(z.poly, x.poly, a)
  return z
end

function mul!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  return mul!(zero(x), x, y)
end

function mul!(
  z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}
) where {T}
  if z === x || z === y
    p = zero(z.poly)
    _mul_p_p!(parent(z), p, x.poly, y.poly)
    z.poly = p
  else
    _mul_p_p!(parent(z), z.poly, x.poly, y.poly)
  end

  return z
end

function mul!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}, a::T) where {T<:FieldElem}
  x.poly = mul!(x.poly, y.poly, a)
  return x
end

function neg!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  x.poly = neg!(x.poly, y.poly)
  return x
end

function pow!(z::PBWAlgebraElem, x::PBWAlgebraElem, n::Int)
  if n == 0
    return one!(z)
  end

  set!(z.poly, x.poly)
  tmp = zero(z)
  for _ in 2:n
    mul!(tmp, z, x)
    swap!(z, tmp)
  end
  return z
end

function sub!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  x.poly = sub!(x.poly, y.poly)
  return x
end

function sub!(z::PBWAlgebraElem{T}, x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}) where {T}
  z.poly = sub!(z.poly, x.poly, y.poly)
  return z
end

function submul!(x::PBWAlgebraElem{T}, y::PBWAlgebraElem{T}, a::T) where {T}
  x.poly = submul!(x.poly, y.poly, a)
  return x
end

###############################################################################
#
#   Multiplication implementation
#
###############################################################################

function _linear_index(A::PBWAlgebra, i::Int, j::Int)
  return div((i - 1) * (2 * ngens(A) - i - 2), 2) + j - 1
end

function _mul_p_p!(
  A::PBWAlgebra{T}, z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::MPolyRingElem{T}
) where {T}
  z = zero!(z)
  cf = zero(coefficient_ring(z))
  mx = Vector{Int}(undef, ngens(parent(z)))
  my = Vector{Int}(undef, ngens(parent(z)))

  for i in 1:length(x)
    exponent_vector!(mx, x, i)
    for j in 1:length(y)
      exponent_vector!(my, y, j)
      cf = mul!(cf, coeff(x, i), coeff(y, j))
      _addmul_m_m!(A, z, mx, my, cf)
    end
  end
end

function _mul_p_m!(
  A::PBWAlgebra{T}, z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::Vector{Int}
) where {T}
  z = zero!(z)
  mx = Vector{Int}(undef, ngens(parent(z)))
  for i in 1:length(x)
    exponent_vector!(mx, x, i)
    _addmul_m_m!(A, z, mx, y, coeff(x, i))
  end
end

function _mul_m_p!(
  A::PBWAlgebra{T}, z::MPolyRingElem{T}, x::Vector{Int}, y::MPolyRingElem{T}
) where {T}
  z = zero!(z)
  my = Vector{Int}(undef, ngens(parent(z)))
  for j in 1:length(y)
    exponent_vector!(my, y, j)
    _addmul_m_m!(A, z, x, my, coeff(y, j))
  end
end

function _addmul_m_m!(
  A::PBWAlgebra{T}, z::MPolyRingElem{T}, x::Vector{Int}, y::Vector{Int}, cf::T
) where {T}
  xl = findlast(!iszero, x)
  if isnothing(xl)
    add_monomial!(z, y, cf)
    return z
  end

  yf = findfirst(!iszero, y)
  if isnothing(yf)
    add_monomial!(z, x, cf)
    return z
  end

  # monomials are ordered, no need for exchange relations
  if xl <= yf
    m = copy(x)
    m .+= y
    add_monomial!(z, m, cf)
    return z
  end

  xll = findprev(!iszero, x, xl - 1)
  yff = findnext(!iszero, y, yf + 1)

  # apply exchange relation
  if isnothing(xll) && isnothing(yff)
    _addmul_gens(A, z, xl, x[xl], yf, y[yf], cf)
    return z
  end

  tmp = zero(z)
  _addmul_gens(A, tmp, xl, x[xl], yf, y[yf], cf)

  m = deepcopy(y)
  tmp2 = zero(z)
  if !isnothing(yff)
    m[yf] = 0
    _mul_p_m!(A, tmp2, tmp, m)
  else
    swap!(tmp2, tmp)
  end
  if !isnothing(xll)
    unsafe_copyto!(m, 1, x, 1, ngens(A))
    m[xl] = 0
    _mul_m_p!(A, tmp, m, tmp2)
  else
    swap!(tmp, tmp2)
  end

  return add!(z, tmp)
end

# i > j
function _addmul_gens(
  A::PBWAlgebra{T}, z::MPolyRingElem{T}, i::Int, n::Int, j::Int, m::Int, cf::T
) where {T}
  ind = _linear_index(A, j, i)
  mon = Vector{Int}(undef, ngens(A))
  for k in 1:ngens(A)
    mon[k] = 0
  end

  # quasi-commutative case
  if length(A.mult[ind][1, 1]) == 1
    mon[i], mon[j] = n, m
    k = add_monomial!(z, mon, mul!(coeff(A.mult[ind][1, 1], 1)^(n * m), cf))
    return z
  end

  # current and required multiplication table size
  curSize = size(A.mult[ind], 1)
  reqSize = max(n, m)
  if (curSize >= reqSize)
    if isassigned(A.mult[ind], n, m)
      return addmul!(z, A.mult[ind][n, m], cf)
    end
  else
    newSize = reqSize + pbwAlg_multGrow
    mult = Matrix{MPolyRingElem{T}}(undef, newSize, newSize)
    for k in 1:curSize
      copyto!(mult, (k - 1) * newSize + 1, A.mult[ind], (k - 1) * curSize + 1, curSize)
    end
    A.mult[ind] = mult
  end

  mon[i] = 1
  for k in 2:n
    if !isassigned(A.mult[ind], k, 1)
      A.mult[ind][k, 1] = zero(A.R)
      _mul_m_p!(A, A.mult[ind][k, 1], mon, A.mult[ind][k - 1, 1])
    end
  end

  mon[i], mon[j] = 0, 1
  for k in 2:m
    if !isassigned(A.mult[ind], n, k)
      A.mult[ind][n, k] = zero(A.R)
      _mul_p_m!(A, A.mult[ind][n, k], A.mult[ind][n, k - 1], mon)
    end
  end

  return addmul!(z, A.mult[ind][n, m], cf)
end

###############################################################################
#
#   Conformance 
#
###############################################################################

function ConformanceTests.generate_element(A::PBWAlgebra)
  return PBWAlgebraElem(A, ConformanceTests.generate_element(A.R))
end
