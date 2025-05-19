function coefficient_ring(R::MPolyRing{T}) where {T}
  return R.coefficient_ring::parent_type(T)
end

function coefficient_ring(x::MPolyRingElem)
  return coefficient_ring(parent(x))
end

function gen(R::MPolyRing, i::Int)
  return MPolyRingElem(R, one(coefficient_ring(R)), [i == j ? 1 : 0 for j in 1:ngens(R)])
end

function gens(R::MPolyRing)
  return [gen(R, i) for i in 1:ngens(R)]
end

function parent(x::MPolyRingElem)
  return x.parent
end

###############################################################################
#
#   Type system
#
###############################################################################

function elem_type(::Type{MPolyRing{T}}) where {T}
  return MPolyRingElem{T}
end

function parent_type(::Type{MPolyRingElem{T}}) where {T}
  return MPolyRing{T}
end

###############################################################################
#
#   
#
###############################################################################

function deepcopy_internal(x::MPolyRingElem, dict::IdDict)
  return get!(dict, x) do
    MPolyRingElem(
      x.parent,
      deepcopy_internal(x.coeffs, dict),
      deepcopy_internal(x.exps, dict),
      deepcopy_internal(x.len, dict),
    )
  end
end

function hash(R::MPolyRing, h::UInt)
  b = 0x890a00575501672e % UInt
  h = hash(R.N, h)

  return xor(h, b)
end

function hash(x::MPolyRingElem, h::UInt)
  b = 0xe445bcf570546ffd % UInt
  h = hash(x.parent, h)
  h = hash(x.coeffs, h)
  h = hash(x.exps, h)
  h = hash(x.len, h)

  return xor(h, b)
end

###############################################################################
#
#   Monomial functions
#
###############################################################################

function monomial_cmp(x::MPolyRingElem, i::Int, m::Memory{Int})
  return monomial_cmp(parent(x), pointer(x.exps, (i - 1) * parent(x).N + 1), pointer(m))
end

function monomial_cmp(R::MPolyRing, m::Memory{Int}, n::Memory{Int})
  for i in 1:(R.N)
    @inbounds cmp = m[i] - n[i]
    if cmp != 0
      return cmp
    end
  end

  return 0
end

function monomial_cmp(R::MPolyRing, m::Ptr{Int}, n::Ptr{Int})
  for i in 1:(R.N)
    cmp = unsafe_load(m, i) - unsafe_load(n, i)
    if cmp != 0
      return cmp
    end
  end

  return 0
end

function monomial_set!(x::MPolyRingElem, i::Int, y::MPolyRingElem, j::Int)
  N = parent(x).N
  unsafe_copyto!(x.exps, (i - 1) * N + 1, y.exps, (j - 1) * N + 1, N)
end

function monomial_set!(x::MPolyRingElem, i::Int, m::Memory{Int})
  N = parent(x).N
  unsafe_copyto!(x.exps, (i - 1) * N + 1, m, 1, N)
end

function add_monomial!(x::MPolyRingElem{T}, m::Memory{Int}) where {T}
  i = 1
  while i <= length(x)
    cmp = monomial_cmp(x, i, m)
    if cmp == 0
      return i
    elseif cmp < 0
      break
    end
    i += 1
  end

  N = parent(x).N
  if length(x.coeffs) > length(x)
    if isnothing(x.coeffs[length(x) + 1])
      x.coeffs[length(x) + 1] = zero(coefficient_ring(x))
    end
    c = x.coeffs[length(x) + 1]

    unsafe_copyto!(x.coeffs, i + 1, x.coeffs, i, length(x) + 1 - i)
    unsafe_copyto!(x.exps, i * N + 1, x.exps, (i - 1) * N + 1, (length(x) + 1 - i) * N)

    x.coeffs[i] = zero!(c)
    unsafe_copyto!(x.exps, (i - 1) * N + 1, m, 1, N)
  else
    coeffs, exps = fit(x, length(x) + 1)

    unsafe_copyto!(coeffs, 1, x.coeffs, 1, i - 1)
    unsafe_copyto!(exps, 1, x.exps, 1, (i - 1) * N)

    coeffs[i] = zero(coefficient_ring(x))
    unsafe_copyto!(exps, (i - 1) * N + 1, m, 1, N)

    unsafe_copyto!(coeffs, i + 1, x.coeffs, i, length(x.coeffs) + 1 - i)
    unsafe_copyto!(exps, i * N + 1, x.exps, (i - 1) * N + 1, (length(x) + 1 - i) * N)

    x.coeffs = coeffs
    x.exps = exps
  end

  x.len += 1
  return i
end

###############################################################################
#
#   Unsafe methods
#
###############################################################################

function fit(x::MPolyRingElem{T}, len::Int) where {T}
  R = parent(x)
  if length(x.coeffs) >= len
    return x.coeffs, x.exps
  end

  if 2 * length(x.coeffs) > len
    len = 2 * length(x.coeffs)
    if length(x.coeffs) >= 256
      cap = length(x.coeffs)
      while true
        cap += (cap + 3 * 256) >> UInt(2)
        if UInt(cap) >= UInt(len)
          break
        end
      end
      if cap >= 0
        len = cap
      end
    end
  end

  return fill!(Memory{Union{T,Nothing}}(undef, len), nothing), Memory{Int}(undef, R.N * len)
end

function fit!(x::MPolyRingElem{T}, len::Int) where {T}
  coeffs, exps = fit(x, len)
  if coeffs !== x.coeffs
    unsafe_copyto!(coeffs, 1, x.coeffs, 1, length(x.coeffs))
    unsafe_copyto!(exps, 1, x.exps, 1, length(x.exps))
    x.coeffs = coeffs
    x.exps = exps
  end
end

function add!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  R = parent(x)

  len = length(x) + length(y)
  fit!(x, len)

  for i in length(x):-1:1
    x.coeffs[length(y) + i] = x.coeffs[i]
    monomial_set!(x, length(y) + i, x, i)
  end

  i = length(y) + 1
  j = k = 1

  while i <= len && j <= length(y)
    cmp = monomial_cmp(
      R, pointer(x.exps, (i - 1) * R.N + 1), pointer(y.exps, (j - 1) * R.N + 1)
    )
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      x.coeffs[i] = nothing
      monomial_set!(x, k, x, i)
      i += 1
    elseif cmp == 0
      x.coeffs[k] = add!(x.coeffs[i], y.coeffs[j])
      if !iszero(x.coeffs[k])
        x.coeffs[i] = nothing
        monomial_set!(x, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      x.coeffs[k] = deepcopy(y.coeffs[j])
      monomial_set!(x, k, y, j)
      j += 1
    end
    k += 1
  end

  while i <= len
    x.coeffs[k] = x.coeffs[i]
    x.coeffs[i] = nothing
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end

  while j <= length(y)
    x.coeffs[k] = deepcopy(y.coeffs[j])
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  return x
end

function addmul!(x::MPolyRingElem{T}, y::MPolyRingElem{T}, a::T) where {T}
  R = parent(x)

  len = length(x) + length(y)
  fit!(x, len)

  for i in length(x):-1:1
    x.coeffs[length(y) + i] = x.coeffs[i]
    monomial_set!(x, length(y) + i, x, i)
  end

  i = length(y) + 1
  j = k = 1

  while i <= len && j <= length(y)
    cmp = monomial_cmp(
      R, pointer(x.exps, (i - 1) * R.N + 1), pointer(y.exps, (j - 1) * R.N + 1)
    )
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      x.coeffs[i] = nothing
      monomial_set!(x, k, x, i)
      i += 1
    elseif cmp == 0
      x.coeffs[k] = addmul!(x.coeffs[i], y.coeffs[j], a)
      if !iszero(x.coeffs[k])
        x.coeffs[i] = nothing
        monomial_set!(x, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      x.coeffs[k] = y.coeffs[j] * a
      monomial_set!(x, k, y, j)
      j += 1
    end
    k += 1
  end

  while i <= len
    x.coeffs[k] = x.coeffs[i]
    x.coeffs[i] = nothing
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end

  while j <= length(y)
    x.coeffs[k] = y.coeffs[j] * a
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  return x
end

function add!(z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  R = parent(z)
  fit!(z, length(x) + length(y))

  i = j = k = 1
  while i <= length(x) && j <= length(y)
    cmp = monomial_cmp(
      R, pointer(x.exps, (i - 1) * R.N + 1), pointer(y.exps, (j - 1) * R.N + 1)
    )
    if cmp > 0
      if isnothing(z.coeffs[k])
        z.coeffs[k] = deepcopy(x.coeffs[i])
      else
        z.coeffs[k] = set!(z.coeffs[k], x.coeffs[i])
      end
      monomial_set!(z, k, x, i)
      i += 1
    elseif cmp == 0
      if isnothing(z.coeffs[k])
        z.coeffs[k] = x.coeffs[i] + y.coeffs[j]
      else
        z.coeffs[k] = add!(z.coeffs[k], x.coeffs[i], y.coeffs[j])
      end
      if !iszero(x.coeffs[k])
        monomial_set!(z, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      if isnothing(z.coeffs[k])
        z.coeffs[k] = deepcopy(y.coeffs[j])
      else
        z.coeffs[k] = set!(z.coeffs[k], y.coeffs[j])
      end
      monomial_set!(z, k, y, j)
      j += 1
    end
    k += 1
  end

  while i <= length(x)
    if isnothing(z.coeffs[k])
      z.coeffs[k] = deepcopy(x.coeffs[i])
    else
      z.coeffs[k] = set!(z.coeffs[k], x.coeffs[i])
    end
    monomial_set!(z, k, x, i)
    i += 1
    k += 1
  end

  while j <= length(y)
    if isnothing(z.coeffs[k])
      z.coeffs[k] = deepcopy(y.coeffs[j])
    else
      z.coeffs[k] = set!(z.coeffs[k], y.coeffs[j])
    end
    monomial_set!(z, k, y, j)
    j += 1
    k += 1
  end

  z.len = k - 1
  return z
end

function mul!(z::MPolyRingElem{T}, x::MPolyRingElem{T}, a::T) where {T}
  if iszero(a)
    return zero!(z)
  elseif isone(a)
    return set!(z, x)
  end

  fit!(z, length(x))
  for i in 1:length(x)
    if isnothing(z.coeffs[i])
      z.coeffs[i] = x.coeffs[i] * a
    else
      z.coeffs[i] = mul!(z.coeffs[i], x.coeffs[i], a)
    end
    monomial_set!(z, i, x, i)
  end

  z.len = length(x)
  return z
end

function set!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  if x === y
    return x
  end

  fit!(x, y.len)
  for i in 1:length(y.coeffs)
    if isnothing(x.coeffs[i])
      x.coeffs[i] = deepcopy(y.coeffs[i])
    else
      x.coeffs[i] = set!(x.coeffs[i], y.coeffs[i])
    end
    monomial_set!(x, i, y, i)
  end

  x.len = y.len
  return x
end

function swap!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  x.coeffs, y.coeffs = y.coeffs, x.coeffs
  x.exps, y.exps = y.exps, x.exps
  x.len, y.len = y.len, x.len
end

###############################################################################
#
#   Arithemtic
#
###############################################################################

function (R::MPolyRing)()
  return zero(R)
end

function (R::MPolyRing{T})(coeff::T, exp::Vector{Int}) where {T}
  return MPolyRingElem(R, coeff, exp)
end

function ngens(R::MPolyRing)
  return R.N
end

function nvars(R::MPolyRing)
  return R.N
end

function one(R::MPolyRing)
  x = zero(R)
  fit!(x, 1)
  x.coeffs[1] = one(coefficient_ring(R))
  return x
end

function iszero(x::MPolyRingElem)
  return x.len == 0
end

function zero(R::MPolyRing)
  return MPolyRingElem(R)
end

function zero!(x::MPolyRingElem)
  x.len = 0
  return x
end

function coeff(x::MPolyRingElem, i::Int)
  return x.coeffs[i]
end

function exponent_vector(x::MPolyRingElem, i::Int)
  R = parent(x)
  m = Vector{Int}(undef, ngens(R))
  copyto!(m, 1, x.exps, (i - 1) * R.N + 1, ngens(R))
  return m
end

function exponent_vector!(exp::Memory{Int}, x::MPolyRingElem, i::Int)
  unsafe_copyto!(exp, 1, x.exps, (i - 1) * parent(x).N + 1, ngens(parent(x)))
end

function length(x::MPolyRingElem)
  return x.len
end

function symbols(R::MPolyRing)
  return R.vars
end

###############################################################################
#
#   Arithemtic
#
###############################################################################

function Base.:+(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  check_parent(x, y)
  return add!(zero(x), x, y)
end
