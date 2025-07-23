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

function is_gen(x::MPolyRingElem)
  fl, _ = is_gen_with_index(x)
  return fl
end

function is_gen_with_index(x::MPolyRingElem)
  if length(x) != 1 || !isone(coeff(x, 1))
    return false, 0
  end

  i = 0
  for j in 1:ngens(parent(x))
    if x.exps[j] != 0
      if i != 0
        return false, 0
      end
      i = j
    end
  end
  return i > 0, i
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

function base_ring(R::MPolyRing{T}) where {T}
  return R.coefficient_ring
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

function monomial_cmp(R::MPolyRing, x::Vector{Int}, i::Int, y::Vector{Int}, j::Int)
  xk = (i - 1) * R.N
  yk = (j - 1) * R.N
  for _ in 1:(R.N)
    xk += 1
    yk += 1
    @inbounds cmp = x[xk] - y[yk]
    if cmp != 0
      return cmp
    end
  end

  return 0
end

function monomial_set!(dst::Vector{Int}, i::Int, src::Vector{Int}, j::Int, N)
  unsafe_copyto!(dst, (i - 1) * N + 1, src, (j - 1) * N + 1, N)
end

function monomial_set!(x::MPolyRingElem, i::Int, y::MPolyRingElem, j::Int)
  monomial_set!(x.exps, i, y.exps, j, parent(x).N)
end

function monomial_set!(x::MPolyRingElem, i::Int, m::Vector{Int})
  monomial_set!(x.exps, i, m, 1, parent(x).N)
end

function add_monomial!(x::MPolyRingElem{T}, m::Vector{Int}, cf::T) where {T}
  i = 1
  while i <= length(x)
    cmp = monomial_cmp(parent(x), x.exps, i, m, 1)
    if cmp == 0
      x.coeffs[i] = add!(x.coeffs[i]::T, cf)
      return nothing
    elseif cmp < 0
      break
    end
    i += 1
  end

  N = parent(x).N
  fit!(x, length(x) + 1)
  if !isnothing(x.coeffs[length(x) + 1])
    cf = set!(x.coeffs[length(x) + 1]::T, cf)
  else
    cf = deepcopy(cf)
  end

  unsafe_copyto!(x.coeffs, i + 1, x.coeffs, i, length(x) + 1 - i)
  unsafe_copyto!(x.exps, i * N + 1, x.exps, (i - 1) * N + 1, (length(x) + 1 - i) * N)

  x.coeffs[i] = cf
  monomial_set!(x, i, m)

  x.len += 1
  return nothing
end

###############################################################################
#
#   Unsafe methods
#
###############################################################################

function fit!(x::MPolyRingElem{T}, len::Int) where {T}
  ol = length(x.coeffs)
  x.coeffs = resize!(x.coeffs, len)
  x.exps = resize!(x.exps, parent(x).N * len)
  for i in (ol + 1):len
    x.coeffs[i] = nothing
  end
end

function add!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  if iszero(x)
    return set!(x, y)
  elseif iszero(y)
    return x
  end
  
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
    cmp = monomial_cmp(R, x.exps, i, y.exps, j)
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      x.coeffs[i] = nothing
      monomial_set!(x, k, x, i)
      i += 1
    elseif cmp == 0
      x.coeffs[k] = add!(x.coeffs[i]::T, y.coeffs[j]::T)
      if !iszero(x.coeffs[k]::T)
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

  if k == i
    x.len = len
    return x
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

function add!(z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  if z === x
    return add!(z, y)
  elseif z === y
    return add!(z, x)
  end
  
  R = parent(z)
  fit!(z, length(x) + length(y))

  i = j = k = 1
  while i <= length(x) && j <= length(y)
    cmp = monomial_cmp(R, x.exps, i, y.exps, j)
    if cmp > 0
      if isnothing(z.coeffs[k])
        z.coeffs[k] = deepcopy(x.coeffs[i])
      else
        z.coeffs[k] = set!(z.coeffs[k]::T, x.coeffs[i]::T)
      end
      monomial_set!(z, k, x, i)
      i += 1
    elseif cmp == 0
      if isnothing(z.coeffs[k])
        z.coeffs[k] = x.coeffs[i]::T + y.coeffs[j]::T
      else
        z.coeffs[k] = add!(z.coeffs[k]::T, x.coeffs[i]::T, y.coeffs[j]::T)
      end
      if !iszero(z.coeffs[k]::T)
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
        z.coeffs[k] = set!(z.coeffs[k]::T, y.coeffs[j]::T)
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
      z.coeffs[k] = set!(z.coeffs[k]::T, y.coeffs[j]::T)
    end
    monomial_set!(z, k, y, j)
    j += 1
    k += 1
  end

  z.len = k - 1
  return z
end

function addmul!(x::MPolyRingElem{T}, y::MPolyRingElem{T}, a::T) where {T}
  if iszero(y) || iszero(a)
    return x
  elseif isone(a)
    return add!(x, y)
  end

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
    cmp = monomial_cmp(R, x.exps, i, y.exps, j)
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      x.coeffs[i] = nothing
      monomial_set!(x, k, x, i)
      i += 1
    elseif cmp == 0
      x.coeffs[k] = addmul!(x.coeffs[i]::T, y.coeffs[j]::T, a)
      if !iszero(x.coeffs[k]::T)
        x.coeffs[i] = nothing
        monomial_set!(x, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      x.coeffs[k] = y.coeffs[j]::T * a
      monomial_set!(x, k, y, j)
      j += 1
    end
    k += 1
  end

  if k == i
    x.len = len
    return x
  end

  while i <= len
    x.coeffs[k] = x.coeffs[i]
    x.coeffs[i] = nothing
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end

  while j <= length(y)
    x.coeffs[k] = y.coeffs[j]::T * a
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  return x
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
      z.coeffs[i] = x.coeffs[i]::T * a
    else
      z.coeffs[i] = mul!(z.coeffs[i]::T, x.coeffs[i]::T, a)
    end
    monomial_set!(z, i, x, i)
  end

  z.len = length(x)
  return z
end

function sub!(z::MPolyRingElem{T}, x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  if z === x
    return sub!(z, y)
  elseif z === y
    return add!(neg!(z), x)
  end
  
  R = parent(z)
  fit!(z, length(x) + length(y))

  i = j = k = 1
  while i <= length(x) && j <= length(y)
    cmp = monomial_cmp(R, x.exps, i, y.exps, j)
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
        z.coeffs[k] = x.coeffs[i] - y.coeffs[j]
      else
        z.coeffs[k] = sub!(z.coeffs[k], x.coeffs[i], y.coeffs[j])
      end
      if !iszero(z.coeffs[k])
        monomial_set!(z, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      if isnothing(z.coeffs[k])
        z.coeffs[k] = -y.coeffs[j]
      else
        z.coeffs[k] = neg!(z.coeffs[k], y.coeffs[j])
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
      z.coeffs[k] = -y.coeffs[j]
    else
      z.coeffs[k] = neg!(z.coeffs[k], y.coeffs[j])
    end
    monomial_set!(z, k, y, j)
    j += 1
    k += 1
  end

  z.len = k - 1
  return z
end

function sub!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
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
    cmp = monomial_cmp(R, x.exps, i, y.exps, j)
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      x.coeffs[i] = nothing
      monomial_set!(x, k, x, i)
      i += 1
    elseif cmp == 0
      x.coeffs[k] = sub!(x.coeffs[i]::T, y.coeffs[j]::T)
      if !iszero(x.coeffs[k]::T)
        x.coeffs[i] = nothing
        monomial_set!(x, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      x.coeffs[k] = -y.coeffs[j]
      monomial_set!(x, k, y, j)
      j += 1
    end
    k += 1
  end

  # if k == i all monomials of y were smaller than those of x
  # since a monomials of x were pushed to the end
  # we are done
  if k == i
    x.len = len
    return x
  end

  while i <= len
    x.coeffs[k] = x.coeffs[i]
    x.coeffs[i] = nothing
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end

  while j <= length(y)
    x.coeffs[k] = -y.coeffs[j]
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  return x
end

function submul!(x::MPolyRingElem{T}, y::MPolyRingElem{T}, a::T) where {T}
  if iszero(y) || iszero(a)
    return x
  elseif isone(a)
    return sub!(x, y)
  end

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
    cmp = monomial_cmp(R, x.exps, i, y.exps, j)
    if cmp > 0
      x.coeffs[k] = x.coeffs[i]
      x.coeffs[i] = nothing
      monomial_set!(x, k, x, i)
      i += 1
    elseif cmp == 0
      x.coeffs[k] = submul!(x.coeffs[i]::T, y.coeffs[j]::T, a)
      if !iszero(x.coeffs[k]::T)
        x.coeffs[i] = nothing
        monomial_set!(x, k, x, i)
      else
        k -= 1
      end

      i += 1
      j += 1
    else
      x.coeffs[k] = submul!(zero(a), y.coeffs[j]::T, a)
      monomial_set!(x, k, y, j)
      j += 1
    end
    k += 1
  end

  if k == i
    x.len = len
    return x
  end

  while i <= len
    x.coeffs[k] = x.coeffs[i]
    x.coeffs[i] = nothing
    monomial_set!(x, k, x, i)
    i += 1
    k += 1
  end

  while j <= length(y)
    x.coeffs[k] = submul!(zero(a), y.coeffs[j]::T, a)
    monomial_set!(x, k, y, j)
    j += 1
    k += 1
  end

  x.len = k - 1
  return x
end

function set!(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  if x === y
    return x
  end

  fit!(x, length(y))
  for i in 1:length(y)
    if isnothing(x.coeffs[i])
      x.coeffs[i] = deepcopy(y.coeffs[i]::T)
    else
      x.coeffs[i] = set!(x.coeffs[i]::T, y.coeffs[i]::T)
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
#   
#
###############################################################################

function (R::MPolyRing)()
  return zero(R)
end

function (R::MPolyRing)(a::Integer)
  if iszero(a)
    return zero(R)
  end

  z = one(R)
  z.coeffs[1] = coefficient_ring(R)(a)
  return z
end

function (R::MPolyRing{T})(a::T) where {T}
  if iszero(a)
    return zero(R)
  end

  z = one(R)
  z.coeffs[1] = a
  return z
end

function (R::MPolyRing{T})(x::MPolyRingElem{T}) where {T}
  @req R === parent(x) "parent mismatch"
  return x
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
  return one!(zero(R))
end

function one!(x::MPolyRingElem{T}) where {T}
  fit!(x, 1)
  if isnothing(x.coeffs[1])
    x.coeffs[1] = one(coefficient_ring(x))
  else
    x.coeffs[1] = one!(x.coeffs[1]::T)
  end
  for i in 1:(parent(x).N)
    x.exps[i] = 0
  end
  x.len = 1
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

function neg!(z::MPolyRingElem{T}, x::MPolyRingElem{T}) where {T}
  fit!(z, length(x))
  for i in 1:length(x)
    if isnothing(z.coeffs[i])
      z.coeffs[i] = -x.coeffs[i]::T
    else
      z.coeffs[i] = neg!(z.coeffs[i]::T, x.coeffs[i]::T)
    end
    monomial_set!(z, i, x, i)
  end

  z.len = length(x)
  return z
end

function coeff(x::MPolyRingElem{T}, i::Int) where {T}
  return x.coeffs[i]::T
end

function setcoeff!(x::MPolyRingElem{T}, i::Int, a::T) where {T}
  x.coeffs[i] = a
  return x
end

function exponent(x::MPolyRingElem, i::Int, j::Int)
  return x.exps[(i - 1) * parent(x).N + j]
end

function set_exponent!(x::MPolyRingElem, i::Int, j::Int, e::Int)
  x.exps[(i - 1) * parent(x).N + j] = e
end

function exponent_vector(x::MPolyRingElem, i::Int)
  return exponent_vector!(Vector{Int}(undef, ngens(parent(x))), x, i)
end

function exponent_vector!(exp::Vector{Int}, x::MPolyRingElem, i::Int)
  R = parent(x)
  unsafe_copyto!(exp, 1, x.exps, (i - 1) * R.N + 1, ngens(R))
  return exp
end

function set_exponent_vector!(x::MPolyRingElem, i::Int, exp::Vector{Int})
  R = parent(x)
  unsafe_copyto!(x.exps, (i - 1) * R.N + 1, exp, 1, R.N)
  return x
end

function length(x::MPolyRingElem)
  return x.len
end

function symbols(R::MPolyRing)
  return R.vars
end

function Base.:(==)(x::MPolyRingElem{T}, y::MPolyRingElem{T}) where {T}
  check_parent(x, y)
  if length(x) != length(y)
    return false
  end

  for i in 1:length(x)
    if monomial_cmp(parent(x), x.exps, i, y.exps, i) != 0
      return false
    end
    if x.coeffs[i] != y.coeffs[i]
      return false
    end
  end

  return true
end

function leading_monomial(x::MPolyRingElem)
  @req !iszero(x) "zero element does not have a leading monomial"
  
  R = parent(x)
  return MPolyRingElem(R, one(coefficient_ring(R)), exponent_vector(x, 1))
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

function divexact(x::MPolyRingElem{T}, a::T) where {T}
  return div!(zero(x), x, a)
end

function div!(z::MPolyRingElem{T}, x::MPolyRingElem{T}, a::T) where {T}
  if iszero(a)
    throw(DivideError())
  elseif isone(a)
    return set!(z, x)
  end

  fit!(z, length(x))
  for i in 1:length(x)
    if isnothing(z.coeffs[i])
      z.coeffs[i] = div(x.coeffs[i]::T, a)
    else
      z.coeffs[i] = div!(z.coeffs[i]::T, x.coeffs[i]::T, a)
    end
    monomial_set!(z, i, x, i)
  end

  z.len = length(x)
  return z
end

###############################################################################
#
#   Conformance 
#
###############################################################################

function (R::MPolyRing{T})(coeffs::Vector{T}, exps::Vector{Vector{Int}}) where {T}
  @req length(coeffs) == length(exps) "lengths must match"
  z = zero(R)
  fit!(z, length(coeffs))
  for i in 1:length(coeffs)
    if iszero(coeffs[i])
      continue
    end
    add_monomial!(z, exps[i], coeffs[i])
  end
  return z
end

function ConformanceTests.generate_element(R::MPolyRing)
  len = rand(0:3)
  coeffs = [ConformanceTests.generate_element(coefficient_ring(R)) for _ in 1:len]
  exps = [[rand(0:3) for _ in 1:ngens(R)] for _ in 1:len]
  return R(coeffs, exps)
end
