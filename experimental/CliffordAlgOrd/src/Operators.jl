
################################################################################
#
#  Unary operators
#
################################################################################

Base.:-(x::Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}) = neg!(parent(x)(), x)

function inv(x::Union{CliffordOrderElem, ZZCliffordOrderElem})
  CO = parent(x)
  CA = ambient_algebra(CO)
  xinv = inv(CA(x))
  if !(xinv in CO)
    error("Element is not invertible")
  end
  return CO(xinv)
end

################################################################################
#
#  Binary operators
#
################################################################################

### Generic ###
function Base.:+(x::T, y::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  check_parent(x, y)
  return parent(x)(coefficients(x) .+ coefficients(y))
end

Base.:-(x::T, y::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}} = x + -y

function Base.:*(x::T, y::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  check_parent(x, y)
  xcoeffs, ycoeffs = copy(coefficients(x)), copy(coefficients(y))
  return parent(x)(_mul_aux(xcoeffs, ycoeffs, gram_matrix(parent(x)), 1))
end

@doc raw"""
    divexact(x::CliffordOrderElem, a::T) where {T<:RingElement} -> CliffordOrderElem
    divexact(x::ZZCliffordOrderElem, a::T) where {T<:RingElement} -> ZZCliffordOrderElem

Return the element `y` in the Clifford order containing `x` such that `ay = x`,
if it exists. If not, an error is raised.
"""
function divexact(x::Union{CliffordOrderElem, ZZCliffordOrderElem}, elt::T) where {T<:RingElement}
  ambalg = algebra(parent(x))
  res = divexact(ambalg(x), elt)
  @req res in parent(x) "Not an exact division"
  return parent(x)(res)
end

##### Algebra #####
@doc raw"""
    divexact(x::CliffordAlgebraElem, a::T) where {T<:RingElement} -> CliffordAlgebraElem

Return the element `y` in the Clifford algebra containing `x` such that `ay = x`,
if it exists. Otherwise an error is raised.
"""
divexact(x::CliffordAlgebraElem, a::T) where {T<:RingElement} =
  parent(x)(divexact.(coefficients(x), a))

##### Order #####
Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U,V}}, a::W) where {T, U<:NumFieldElem, V, W<:FieldElem} = parent(x)(a * algebra(parent(x))(x))
Base.:*(a::W, x::CliffordOrderElem{T, CliffordAlgebra{U,V}}) where {T, U<:NumFieldElem, V, W<:FieldElem} = parent(x)(a * algebra(parent(x))(x))
Base.:*(x::CliffordOrderElem{T, CliffordAlgebra{U,V}}, a::Rational{Int}) where {T, U<:NumFieldElem, V} = parent(x)(a .* coefficients(x))
Base.:*(a::Rational{Int}, x::CliffordOrderElem{T, CliffordAlgebra{U,V}}) where {T, U<:NumFieldElem, V} = parent(x)(a .* coefficients(x))

### ZZ ###
Base.:*(x::ZZCliffordOrderElem, a::QQFieldElem) = parent(x)(a .* coefficients(x))
Base.:*(a::QQFieldElem, x::ZZCliffordOrderElem) = parent(x)(a .* coefficients(x))
Base.:*(x::ZZCliffordOrderElem, a::Rational{Int}) = parent(x)(a .* coefficients(x))
Base.:*(a::Rational{Int}, x::ZZCliffordOrderElem) = parent(x)(a .* coefficients(x))

################################################################################
#
#  Auxillary function for multiplication 
#
################################################################################

function _mul_aux(x::Vector{T}, y::Vector{T}, gram::MatElem{T}, i::Int) where {T<:RingElement}
  R = base_ring(gram)
  out = [zero(R) for _ in 1:length(x)]
  temp_buffers = [[zero(R) for _ in 1:length(x)] for _ in 1:ncols(gram)]
  _mul_aux!(out, x, y, gram, i, temp_buffers)
  return out
end

function _mul_aux!(out::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}, gram::MatElem{T}, i::Int, temp_buffers::Vector{Vector{T}}) where {T<:RingElement}
  if is_zero(y)
    return out
  elseif length(y) == 1
    y_val = y[1]
    if !is_zero(y_val)
      @inbounds for j in eachindex(out, x)
        addmul!(out[j], x[j], y_val)
      end
    end
    return out
  end
    
  y_even = @view y[1:2:end]
  y_odd  = @view y[2:2:end]
    
  # Accumulate x * y_even
  _mul_aux!(out, x, y_even, gram, i + 1, temp_buffers)
    
  # Accumulate (x * e_i) * y_odd
  if !is_zero(y_odd)
    x_next = temp_buffers[i]
    _mul_with_gen!(x_next, x, i, gram)
    _mul_aux!(out, x_next, y_odd, gram, i + 1, temp_buffers)
  end
    
  return out
end

# Implements right multiplication with the 'i'-th generator.
function _mul_with_gen!(out::AbstractVector{T}, x::AbstractVector{T}, i::Int, gram::MatElem{T}) where {T<:RingElement}
  # Reset the buffer
  @inbounds for j in eachindex(out)
    zero!(out[j])
  end
  
  @inbounds for char in 1:length(x)
    if !is_zero(x[char])
      _add_mul_baseelt_with_gen!(out, x[char], char, i, gram, 0)
    end
  end
  return out
end

function _add_mul_baseelt_with_gen!(out::AbstractVector{T}, c::T, char::Int, i::Int, gram::MatElem{T}, shift::Int=0) where {T <: RingElement}
  is_positive = true

  while true
    if char == 1
      idx = 1 + (1 << (i - 1)) + shift
      is_positive ? add!(out[idx], out[idx], c) : sub!(out[idx], out[idx], c)
      return
    end
        
    j = 8 * sizeof(Int) - leading_zeros(char - 1)
        
    if j < i
      idx = char + (1 << (i - 1)) + shift
      is_positive ? add!(out[idx], out[idx], c) : sub!(out[idx], out[idx], c)
      return
    end
        
    if j == i
      idx = char - (1 << (i - 1)) + shift
      half_G = divexact(gram[i, i], 2)
      is_positive ? addmul!(out[idx], c, half_G) : submul!(out[idx], c, half_G)
      return
    end
        
    # j > i
    idx = char - (1 << (j - 1)) + shift
    is_positive ? addmul!(out[idx], c, gram[i, j]) : submul!(out[idx], c, gram[i, j])
    
    # Setup for the next loop iteration
    is_positive = !is_positive
    char = char - (1 << (j - 1))
    shift = shift + (1 << (j - 1))
  end
end

################################################################################
#
#  Equality and hash
#
################################################################################

Base.:(==)(x::T, y::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}} = parent(x) === parent(y) && coefficients(x) == coefficients(y)

function Base.hash(x::CliffordAlgebraElem, h::UInt)
  b = 0x1c4629b4de23b24c % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

function Base.hash(x::CliffordOrderElem, h::UInt)
  b = 0x8f3a7c2b1d4e5f6a % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

function Base.hash(x::ZZCliffordOrderElem, h::UInt)
  b = 0x924e7a492844d7a0 % UInt
  h = hash(parent(x), h)
  h = hash(coefficients(x), h)
  return xor(h, b)
end

################################################################################
#
#  Unsafe operations 
#
################################################################################

function zero!(c::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  cs = coefficients(c)
  for i in 1:length(cs)
    zero!(cs[i])
  end
  return c
end

function one!(c::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  cs = coefficients(c)
  one!(cs[1])
  for i in 2:length(cs)
    zero!(cs[i])
  end
  return c
end

function neg!(c::T, a::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  cs, as = coefficients(c), coefficients(a)
  for i in 1:length(cs)
    neg!(cs[i], as[i])
  end 
  return c
end

function add!(c::T, a::T, b::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  as, bs, cs = coefficients(a), coefficients(b), coefficients(c)
  for i in 1:length(cs)
    add!(cs[i], as[i], bs[i])
  end
  return c
end

function sub!(c::T, a::T, b::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  as, bs, cs = coefficients(a), coefficients(b), coefficients(c)
  for i in 1:length(cs)
    sub!(cs[i], as[i], bs[i])
  end
  return c
end

mul!(a::T, b::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}} = mul!(a, a, b)

function mul!(c::T, a::T, b::T) where {T <: Union{CliffordAlgebraElem, CliffordOrderElem, ZZCliffordOrderElem}}
  cs, as, bs = coefficients(c), coefficients(a), coefficients(b)
  
  if cs === as || cs === bs
    tmp = a * b
    tmps = coefficients(tmp)
    @inbounds for i in eachindex(cs)
      zero!(cs[i])
      add!(cs[i], cs[i], tmps[i])
    end
    return c
  end

  gram = gram_matrix(parent(c))
  R = base_ring(gram)
    
  temp_buffers = [[R() for _ in 1:length(cs)] for _ in 1:ncols(gram)]
    
  @inbounds for i in eachindex(cs)
    zero!(cs[i])
  end
    
  _mul_aux!(cs, as, bs, gram, 1, temp_buffers)
    
  return c
end
