
################################################################################
#
#  Construction
#
################################################################################

### Algebra ###

@doc raw"""
    clifford_algebra(qs::QuadSpace) -> CliffordAlgebra

Return the Clifford algebra of the quadratic space `qs`.

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, QQ[0 1; 1 0]))
Clifford algebra of quadratic space with Gram matrix
  [0   1]
  [1   0]
defined over rational field

julia> K, a = quadratic_field(-5); C = clifford_algebra(quadratic_space(K, K[2 a; a 2]))
Clifford algebra of quadratic space with Gram matrix
  [       2   sqrt(-5)]
  [sqrt(-5)          2]
defined over imaginary quadratic field defined by x^2 + 5
```
"""
clifford_algebra(qs::Hecke.QuadSpace) =
  CliffordAlgebra{elem_type(base_ring(qs)),typeof(gram_matrix(qs))}(qs)

### Elements ###
(C::CliffordAlgebra)() = CliffordAlgebraElem(C)

function (C::CliffordAlgebra)(a::R) where {R<:RingElem}
  res = fill(zero(a), dim(C))
  res[1] = a
  return CliffordAlgebraElem(C, res)
end

# for disambiguation
function (C::CliffordAlgebra)(a::ZZRingElem)
  res = fill(zero(a), dim(C))
  res[1] = a
  return CliffordAlgebraElem(C, res)
end

function (C::CliffordAlgebra)(a::CliffordAlgebraElem)
  @req parent(a) === C "The element does not lie in the Clifford algebra"
  return a
end

(C::CliffordAlgebra)(coeff::Vector{R}) where {R} = CliffordAlgebraElem(C, coeff)

################################################################################
#
#  Basic field access
#
################################################################################

### Algebra ###
@doc raw"""
    base_ring(C::CliffordAlgebra) -> Ring

Return the base ring of the Clifford algebra `C`.
"""
base_ring(C::CliffordAlgebra) = C.base_ring::base_ring_type(typeof(C))

@doc raw"""
    space(C::CliffordAlgebra) -> QuadSpace

Return the underlying quadratic space of the Clifford algebra `C`.
"""
space(C::CliffordAlgebra{T, S}) where {T, S} = C.space::Hecke.QuadSpace{parent_type(T), S}

@doc raw"""
    gram_matrix(C::CliffordAlgebra) -> MatElem

Return the Gram matrix of the free quadratic module with Clifford algebra `C`.
"""
gram_matrix(C::CliffordAlgebra) = C.gram

@doc raw"""
    dim(C::CliffordAlgebra) -> Int

Return the dimension of the Clifford algebra `C` over its base ring.
"""
dim(C::CliffordAlgebra) = C.dim

### Elements ###
@doc raw"""
    parent(x::CliffordAlgebraElem) -> CliffordAlgebra

Return the Clifford algebra containing `x`.
"""
parent(x::CliffordAlgebraElem) = x.parent

@doc raw"""
    coefficients(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of `x` wrt the
canonical basis of its parent Clifford algebra.  
"""
coefficients(x::CliffordAlgebraElem) = x.coeffs

################################################################################
#
#  basic functionality
#
################################################################################

base_ring_type(::Type{CliffordAlgebra{T, S}}) where {T, S} = parent_type(T)
elem_type(::Type{CliffordAlgebra{T, S}}) where {T, S} = CliffordAlgebraElem{T, S}
is_domain_type(::Type{CliffordAlgebraElem{T, S}}) where {T, S} = false
is_exact_type(::Type{CliffordAlgebraElem{T, S}}) where {T, S} = true
parent_type(::Type{CliffordAlgebraElem{T, S}}) where {T, S} = CliffordAlgebra{T, S}

@doc raw"""
    zero(C::CliffordAlgebra) -> CliffordAlgebraElem

Return the additive identity of the Clifford algebra `C`.
"""
zero(C::CliffordAlgebra) = C()

@doc raw"""
    one(C::CliffordAlgebra) -> CliffordAlgebraElem

Return the multiplicative identity of the Clifford algebra `C`.
"""
one(C::CliffordAlgebra) = C(1)

@doc raw"""
    basis(C::CliffordAlgebra, i::Int) -> CliffordAlgebraElem

Return the `i`-th canonical basis vector of the Clifford algebra `C`.

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, QQ[0 1; 1 0]));

julia> basis(C, 1)
[1, 0, 0, 0]

julia> basis(C, 2)
[0, 1, 0, 0]

julia> basis(C, 3)
[0, 0, 1, 0]

julia> basis(C, 4)
[0, 0, 0, 1]
```
"""
function basis(C::CliffordAlgebra, i::Int)
  R = base_ring(C)
  coeffs = fill(zero(R), dim(C))
  coeffs[i] = one(R)
  return C(coeffs)
end

@doc raw"""
    gen(C::CliffordAlgebra, i::Int) -> CliffordAlgebraElem

Return the `i`-th canonical algebra generator of the Clifford
algebra `C`. This is just the image of the `i`-th basis vector
of the underlying quadratic space under the canonical embedding
into the Clifford algebra.

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, identity_matrix(QQ,3)));

julia> gen(C, 1)
[0, 1, 0, 0, 0, 0, 0, 0]

julia> gen(C, 2)
[0, 0, 1, 0, 0, 0, 0, 0]

julia> gen(C, 3)
[0, 0, 0, 0, 1, 0, 0, 0]
```
"""
function gen(C::CliffordAlgebra, i::Int)
  R = base_ring(C)
  coeffs = fill(zero(R), dim(C))
  i > 0 || throw(BoundsError(coeffs, i))
  coeffs[2^(i - 1) + 1] = one(R)
  return C(coeffs)
end

@doc raw"""
    basis(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the canonical basis of the Clifford algebra `C`.

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, QQ[0 1; 1 0]));

julia> basis(C)
4-element Vector{CliffordAlgebraElem{QQFieldElem, QQMatrix}}:
 [1, 0, 0, 0]
 [0, 1, 0, 0]
 [0, 0, 1, 0]
 [0, 0, 0, 1]
```
"""
basis(C::CliffordAlgebra) = [basis(C, i) for i in 1:dim(C)]

@doc raw"""
    gens(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return the vector of canonical algebra generators of the Clifford algebra `C`, i.e.,
if `gram_matrix(C)` is the Gram matrix of the underlying quadratic space with respect
to the basis (e_1,...,e_n) then the vector of images under the canonical embedding into
the Clifford algebra is returned.

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, identity_matrix(QQ, 3)));

julia> gens(C)
3-element Vector{CliffordAlgebraElem{QQFieldElem, QQMatrix}}:
 [0, 1, 0, 0, 0, 0, 0, 0]
 [0, 0, 1, 0, 0, 0, 0, 0]
 [0, 0, 0, 0, 1, 0, 0, 0]
```
"""
gens(C::CliffordAlgebra) = [gen(C, i) for i in 1:dim(space(C))]

@doc raw"""
    is_commutative(C::CliffordAlgebra) -> Bool

Return `true` if `C` is commutative and `false` otherwise.
"""
is_commutative(C::CliffordAlgebra) = dim(C) == 1 || dim(C) == 2

################################################################################
#
#  Element Access
#
################################################################################

@doc raw"""
    coeff(x::CliffordAlgebraElem, i::Int) -> FieldElem

Return the `i`-th coefficient of the element `x`.
"""
coeff(x::CliffordAlgebraElem, i::Int) = coefficients(x)[i]

getindex(x::CliffordAlgebraElem, i::Int) = coefficients(x)[i]

setindex!(x::CliffordAlgebraElem, newentry::RingElement, i::Int) = (coefficients(x)[i] = base_ring(x)(newentry))

################################################################################
#
#  Other functionality
#
################################################################################

@doc raw"""
    basis_of_centroid(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return a basis of the centroid of `C`. Unless `dim(space(C)) = 0`, it consists of
two elements. The first one is the multiplicative identity of `C`. The square of
the second basis element, if present, equals `quadratic_discriminant(C)`; see
[`quadratic_discriminant`](@ref quadratic_discriminant(C::CliffordOrder)).

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, gram_matrix(root_lattice(:A, 4))));

julia> basis_of_centroid(C)
2-element Vector{CliffordAlgebraElem{QQFieldElem, QQMatrix}}:
 [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 [1, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 4]
```
"""
function basis_of_centroid(C::CliffordAlgebra)
  dim(space(C)) == 0 && return [one(C)]
  return [one(C), _orth_elt(C)]
end

@doc raw"""
    quadratic_discriminant(C::CliffordAlgebra) -> FieldElem

Return the quadratic discriminant of `C` as an element of `base_ring(C)`.
"""
function quadratic_discriminant(C::CliffordAlgebra)
  !isdefined(C, :disq) && _orth_elt(C)
  return C.disq
end

@doc raw"""
    disq(C::CliffordAlgebra) -> FieldElem

Alias for [`quadratic_discriminant`](@ref quadratic_discriminant(C::CliffordAlgebra)).
"""
disq(C::CliffordAlgebra) = quadratic_discriminant(C)

@doc raw"""
    basis_of_center(C::CliffordAlgebra) -> Vector{CliffordAlgebraElem}

Return a basis of the center of `C`. It equals `basis_of_centroid(C)`, if and only
if `dim(space(C))` is odd. Otherwise it contains only the multiplicative identity of `C`.
"""
function basis_of_center(C::CliffordAlgebra)
  is_odd(dim(space(C))) && return basis_of_centroid(C)
  return [one(C)]
end

@doc raw"""
    representation_matrix(x::CliffordAlgebraElem, action::Symbol = :left) -> MatElem

Return the representation matrix of the element `x` with respect to `basis(parent(x))`. The multiplication is from
the left if `action == :left` and from the right if `action == :right`.
"""
function representation_matrix(x::CliffordAlgebraElem, action::Symbol = :left)
  @req (action == :left) || (action == :right) "The action must be :left or :right."
  C = parent(x)
  n = dim(C)
  res = zero_matrix(base_ring(C), n, n)
  if action == :left
    for i in 1:n
      res[i, :] = coefficients(x * basis(C, i))
    end
  elseif action == :right
    for i in 1:n
      res[i, :] = coefficients(basis(C, i) * x)
    end
  end
  return res
end

################################################################################
#
#  Graded parts
#
################################################################################

@doc raw"""
    even_coefficients(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of `x` with respect to the
canonical basis of its parent Clifford algebra,
but with all its coefficients that correspond to basis
elements with odd grading set to zero.
"""
even_coefficients(x::CliffordAlgebraElem) = [iseven(count_ones(y - 1)) ? coefficients(x)[y] : zero(base_ring(parent(x))) for y in 1:dim(parent(x))]

@doc raw"""
    even_part(x::CliffordAlgebraElem) -> CliffordAlgebraElem

Return the projection of `x` onto the even Clifford algebra

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, identity_matrix(QQ,3))); x = C(collect(1:8))
[1, 2, 3, 4, 5, 6, 7, 8]

julia> even_part(x)
[1, 0, 0, 4, 0, 6, 7, 0]
```
"""
even_part(x::CliffordAlgebraElem) = parent(x)(even_coefficients(x))

@doc raw"""
    is_even(x::CliffordAlgebraElem) -> Bool

Return 'true' if 'x' is even, i.e. if [`even_part(x)`](@ref even_part(x::CliffordAlgebraElem))
and `x` coincide. Otherwise, return `false`.
"""
is_even(x::CliffordAlgebraElem) = (x == even_part(x))

@doc raw"""
    odd_coefficients(x::CliffordAlgebraElem) -> Vector

Return the coefficient vector of `x` with respect to the
canonical basis of its parent Clifford algebra,
but with all its coefficients that correspond to basis elements
with even grading set to zero.
"""
odd_coefficients(x::CliffordAlgebraElem) = [isodd(count_ones(y - 1)) ? coefficients(x)[y] : zero(base_ring(parent(x))) for y in 1:dim(parent(x))]

@doc raw"""
    odd_part(x::CliffordAlgebraElem) -> CliffordAlgebraElem

Return the projection of `x` onto the odd Clifford algebra.

# Examples

```jldoctest
julia> C = clifford_algebra(quadratic_space(QQ, identity_matrix(QQ,3))); x = C(collect(1:8))
[1, 2, 3, 4, 5, 6, 7, 8]

julia> odd_part(x)
[0, 2, 3, 0, 5, 0, 0, 8]
```
"""
odd_part(x::CliffordAlgebraElem) = parent(x)(odd_coefficients(x))

@doc raw"""
    is_odd(x::CliffordAlgebraElem) -> Bool

Return 'true' if 'x' is odd, i.e. if [`odd_part(x)`](@ref odd_part(x::CliffordAlgebraElem))
and `x` coincide. Otherwise, return `false`.
"""
is_odd(x::CliffordAlgebraElem) = (x == odd_part(x))

################################################################################
#
#  Promotion rules
#
################################################################################

function AbstractAlgebra.promote_rule(::Type{CAE}, ::Type{V}) where
        {T<:RingElement, S, CAE<:CliffordAlgebraElem{T, S}, V<:RingElement}
    AbstractAlgebra.promote_rule(T, V) == T ? CAE : Union{}
end

AbstractAlgebra.promote_rule(::Type{CAE}, ::Type{CAE}) where
        {T<:RingElement, S, CAE<:CliffordAlgebraElem{T, S}} = CAE

################################################################################
#
#  Random generation 
#
################################################################################

function rand(rng::Random.AbstractRNG, C::CliffordAlgebra, v...)
  coeffs = [rand(rng, base_ring(C), v...) for _ in 1:C.dim]
  return C(coeffs)
end

rand(C::CliffordAlgebra, v...) = rand(Random.default_rng(), C, v...)

################################################################################
#
#  Conformance test element generation 
#
################################################################################

ConformanceTests.generate_element(C::CliffordAlgebra) = C([base_ring(C)(rand(-3:3)) for _ in 1:dim(C)])

################################################################################
#
#  Auxillary functions
#
################################################################################

function _orth_elt(C::CliffordAlgebra)
  isdefined(C, :orth_elt) && isdefined(C, :disq) && return C.orth_elt::elem_type(C)
    
  n = dim(space(C))
  if n == 0
    C.orth_elt = C(1)
    C.disq = base_ring(C)(1)
    return C.orth_elt::elem_type(C)
  end
    
  T = orthogonal_basis(space(C))  
  orth_gens = [sum(gen(C, j) * T[i, j] for j in 1:n) for i in 1:n]
  orth_elt = prod(orth_gens)
  scale = denominator(orth_elt)
  C.orth_elt = orth_elt * scale
  
  if !isdefined(C, :disq)
    sign = n % 4 in (0, 1) ? 1 : -1
    C.disq = sign * prod(coefficients(x^2)[1] for x in orth_gens) * scale^2
  end
  
  return C.orth_elt::elem_type(C)
end

#=
function _mul_aux(x::Vector{T}, y::Vector{T}, gram::MatElem{T}, i::Int) where {T<:RingElement}
  if length(y) == 1
    return x .* y[1]
  elseif is_zero(y)
    return fill(y[1], 2^ncols(gram))
  end
  return _mul_aux(_mul_with_gen(x, i, gram), y[2:2:end], gram, i + 1) +
         _mul_aux(x, y[1:2:end], gram, i + 1)
end

#Implements the right multiplication with 'i'-th generator e_i of the Clifford algebra containing x.
#The result is returned as a coefficient vector for further computations.
_mul_with_gen(x::Vector{T}, i::Int, gram::MatElem{T}) where {T<:RingElement} = sum(
  map(
    char ->
      x[char] .*
      _mul_baseelt_with_gen(char, i, gram),
      1:2^ncols(gram),
  ),
)

#Shift all entries of the vector 'X' to the right by 's' entries without changing the length of X.
function _shift_entries!(X::Vector, s::Int)
  X[(s + 1) : end] = X[1 : (end - s)]
  X[1 : s] = fill(parent(X[1])(), s)
  return X
end

#Right multiplication of the basis element represented by 'char' with the
#'i'-th generator of the Clifford algebra/order containing said basis vector
function _mul_baseelt_with_gen(char::Int, i::Int, gram::MatElem)
  R = base_ring(gram)
  res = fill(R(), 2^ncols(gram))
  if char == 1
    res[char + 2^(i - 1)] = R(1)
    return res
  end
  j = 8 * sizeof(Int) - leading_zeros(char - 1)
  if j < i
    res[char + 2^(i - 1)] = R(1)
    return res
  end
  if j == i
    res[char - 2^(i - 1)] = R(divexact(gram[i, i], 2))
    return res
  end
  res[char - 2^(j - 1)] = gram[i, j]
  res -= _shift_entries!(_mul_baseelt_with_gen(char - 2^(j - 1), i, gram), 2^(j - 1))
end
=#

# Wrapper: Sets up the pre-allocated buffers exactly ONCE
function _mul_aux(x::Vector{T}, y::Vector{T}, gram::MatElem{T}, i::Int) where {T<:RingElement}
  R = base_ring(gram)
  out = [zero(R) for _ in 1:length(x)]
  temp_buffers = [[zero(R) for _ in 1:length(x)] for _ in 1:ncols(gram)]
  _mul_aux!(out, x, y, gram, i, temp_buffers)
  return out
end

# The allocation-free recursive multiplication engine.
function _mul_aux!(out::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}, gram::MatElem{T}, i::Int, temp_buffers::Vector{Vector{T}}) where {T<:RingElement}
  if is_zero(y)
    return out
  elseif length(y) == 1
    y_val = y[1]
    if !is_zero(y_val)
      @inbounds for j in eachindex(out, x)
        out[j] += x[j] * y_val
      end
    end
    return out
  end
    
  # @views prevent memory allocation when slicing the array
  y_even = @view y[1:2:end]
  y_odd  = @view y[2:2:end]
    
  # Accumulate x * y_even
  _mul_aux!(out, x, y_even, gram, i + 1, temp_buffers)
    
  # Accumulate (x * e_i) * y_odd
  if !is_zero(y_odd)
    x_next = temp_buffers[i] # Grab the pre-allocated buffer for this depth
    _mul_with_gen!(x_next, x, i, gram)
    _mul_aux!(out, x_next, y_odd, gram, i + 1, temp_buffers)
  end
    
  return out
end

# Implements right multiplication with the 'i'-th generator.
function _mul_with_gen!(out::AbstractVector{T}, x::AbstractVector{T}, i::Int, gram::MatElem{T}) where {T<:RingElement}
  # Reset the buffer. Safe from mutation bugs because `+=` rebinds the index.
  fill!(out, zero(base_ring(gram))) 
  @inbounds for char in 1:length(x)
    if !is_zero(x[char])
      _add_mul_baseelt_with_gen!(out, x[char], char, i, gram, 0)
    end
  end
  return out
end

# Resolves basis element multiplication WITHOUT creating any arrays or recursive stack calls.
function _add_mul_baseelt_with_gen!(out::AbstractVector, c, char::Int, i::Int, gram::MatElem, shift::Int=0)
  while true
    if char == 1
      @inbounds out[1 + (1 << (i - 1)) + shift] += c
      return
    end
        
    j = 8 * sizeof(Int) - leading_zeros(char - 1)
    if j < i
      @inbounds out[char + (1 << (i - 1)) + shift] += c
      return
    end
    if j == i
      @inbounds out[char - (1 << (i - 1)) + shift] += c * divexact(gram[i, i], 2)
      return
    end
        
    @inbounds out[char - (1 << (j - 1)) + shift] += c * gram[i, j]        
    # Setup the variables for the next loop iteration instead of recursing
    c = -c
    char = char - (1 << (j - 1))
    shift = shift + (1 << (j - 1))
  end
end
