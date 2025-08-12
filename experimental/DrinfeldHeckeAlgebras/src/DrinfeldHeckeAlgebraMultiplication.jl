#######################################
# Methods to multiply elements in Drinfeld-Hecke algebra
#
# Cassandra Koenen, 2025
#######################################

#######################################
# Multiplication methods to break recursion
#######################################

# Multiply with scalar from left and right
function *(c::S, a::DrinfeldHeckeAlgebraElem{T, S}) where {T <: FieldElem, S <: RingElem}
  A = a.parent

  return A(base_algebra(A)(c) * a.element)
end

*(a::DrinfeldHeckeAlgebraElem{T, S}, c::S) where {T <: FieldElem, S <: RingElem} = c * a

# Multiply group element g from the right
function *(a::DrinfeldHeckeAlgebraElem, g::MatrixGroupElem)
  A = a.parent

  return A(a.element * group_algebra(A)(g))
end

#######################################
# Multiplication methods to implement Drinfeld-Hecke algebra multiplication
#######################################

# Multiply general elements a and b
function multiply(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem)
  A = a.parent
  
  if A != b.parent
   throw(ArgumentError("Elements do not belong to the same algebra"))
  end

  if is_one(a) || is_zero(b) return b end
  if is_one(b) || is_zero(a) return a end
  
  (c, m, g, tail) = split_element(b)

  # a * b = a * (c * m * g + tail) = c * (a * m) * g + a * tail
  return c * multiply_a_with_m(a, m) * g + multiply(a, tail)
end

# Multiply Drinfeld-Hecke algebra element a with monomial m
function multiply_a_with_m(a::DrinfeldHeckeAlgebraElem, m::DrinfeldHeckeAlgebraElem)
  if is_one(m) || is_zero(a) return a end

  (x,mm) = split_monomial_left(m)
  
  # a * m = (a * x) * mm
  return multiply_a_with_m(multiply_a_with_x(a, x), mm)
end

# Multiply Drinfeld-Hecke algebra element a with generator x
function multiply_a_with_x(a::DrinfeldHeckeAlgebraElem, x::DrinfeldHeckeAlgebraElem)
  if is_zero(a) return a end

  A = a.parent
  (c, m, g, tail) = split_element(a)
 
  # Let g act on x
  gx = left_action_on_generator(g,x)
  
  # a * x 
  # = (c * m * g + tail) * x 
  # = c * m * g * x    +   tail * x 
  # = c * m * gx * g   +   tail * x
  return c * multiply_m_with_f(m, A(gx)) * g + multiply_a_with_x(tail, x)
end

# Multiply monomial m with polynomial f
function multiply_m_with_f(m::DrinfeldHeckeAlgebraElem, f::DrinfeldHeckeAlgebraElem)
  if is_one(m) || is_zero(f) return f end

  A = m.parent
  (lc, lm, tail) = split_polynomial(f)
  
  # m * f = m * (lc * lm + tail) = lc * m * lm + m * tail
  return lc * multiply_m1_with_m2(m, lm) + multiply_m_with_f(m, tail)
end

# Multiply monomial m1 with monomial m2
function multiply_m1_with_m2(m1::DrinfeldHeckeAlgebraElem, m2::DrinfeldHeckeAlgebraElem)
  A = m1.parent
  
  # cases of empty monomials
  if is_one(m1) return m2 end
  if is_one(m2) return m1 end

  (mm1,x) = split_monomial_right(m1)
  (y,mm2) = split_monomial_left(m2)
  
  if index_of_generator(x) < index_of_generator(y)
    # In this case x and y are in correct order and we can multiply m1 and m2 in the group algebra
    return A(m1.element * m2.element)
  end
  
  # Otherwise a correctional term will come in when multiplying x and y
  κ = form(A)
  xy = A(y.element * x.element) - κ(y,x)
  (c, m, g, tail) = split_element(xy)

  # m1 * m2 = mm1 * (x * y) * mm2 = mm1 * (c * m * g + tail) * mm2 = c * mm1 * m * g * mm2 + mm1 * tail * mm2
  return c * multiply(multiply_m1_with_m2(mm1, m) * g, mm2) + multiply(mm1, multiply_a_with_m(tail, mm2))
end

#######################################
# Helper functions for recursion
#######################################

# Returns 4-tuple (c, m, g, tail) where
# - c is a scalar
# - m is a monomial
# - g a group element
# - tail = a - c * m * g
function split_element(a::DrinfeldHeckeAlgebraElem)
  A = a.parent
  RG = group_algebra(A)
  elm = a.element
  
  for (i,f) in enumerate(elm.coeffs)
    if is_zero(f)
      continue
    end

    c = leading_coefficient(f)
    m = A(leading_monomial(f))
    g = RG.base_to_group[i]
    tail = a - c * m * g

    return (c, m, g, tail)
  end

  throw(ArgumentError("zero element can not be split"))
end

# Returns triple (lc, lm, tail) where
# - lc is the leading term of f
# - lm is the leading monomial of f
# - tail = f - lc * lm
function split_polynomial(f::DrinfeldHeckeAlgebraElem)
  A = f.parent
  elm = f.element.coeffs[1]
  
  if is_zero(elm)
    throw(ArgumentError("zero polynomial can not be split"))
  end

  lc = leading_coefficient(elm)
  lm = A(leading_monomial(elm))
  t = A(tail(elm))

  return (lc, lm, t)
end

# Returns tuple (x, m) where
# - x is a generator of the underlying polynomial ring
# - a = x * m
function split_monomial_left(a::DrinfeldHeckeAlgebraElem)
  A = a.parent
  m = a.element.coeffs[1]
  
  if !is_monomial(m)
    throw(ArgumentError("Error: element " * print(m) * " is not a monomial"))
  end

  for (i, exp) in enumerate(collect(exponents(m))[1])
    if exp > 0
      x = gen(m.parent, i)
      return (A(x), A(m / x))
    end
  end

  throw(ArgumentError("Error: Can not split empty monomial"))
end

# Returns tuple (m, x) where
# - x is a generator of the underlying polynomial ring
# - a = m * x
function split_monomial_right(a::DrinfeldHeckeAlgebraElem)
  A = a.parent
  m = a.element.coeffs[1]
  
  if !is_monomial(m)
    throw(ArgumentError("Error: element " * print(m) * " is not a monomial"))
  end

  reversed_exponents = reverse(collect(exponents(m))[1])
  reversed_generators = reverse(gens(m.parent))
  
  for (i, exp) in enumerate(reversed_exponents)
    if exp > 0
      x = reversed_generators[i]
      return (A(m / x), A(x))
    end
  end

  throw(ArgumentError("Error: Can not split empty monomial"))
end

#######################################
# Helper functions for multiplication
#######################################

(κ::DrinfeldHeckeForm)(x::DrinfeldHeckeAlgebraElem, y::DrinfeldHeckeAlgebraElem) = x.parent(κ(generator_to_vector(x), generator_to_vector(y)))

function generator_to_vector(x::DrinfeldHeckeAlgebraElem)
  A = x.parent
  R = base_ring(A)
  n = degree(group(A))
  v = [R(0) for _ in 1:n]
  v[index_of_generator(x)] = R(1)
  
  return v
end

function vector_to_algebra_element(A::DrinfeldHeckeAlgebra{T, S}, v::Vector{S}) where {T <: FieldElem, S <: RingElem}
  RV = base_algebra(A)
  res = RV()
  
  for i in 1:length(v)
    res = res + v[i] * RV[i]
  end
  
  return A(res)
end

function index_of_generator(x::DrinfeldHeckeAlgebraElem)
  A = x.parent
  n = degree(group(A))
  
  for i in 1:n
    if A(gen(base_algebra(A), i)) == x
      return i
    end
  end

  throw(ArgumentError("Element is not a generator."))
end

function generator_to_polynomial(x::DrinfeldHeckeAlgebraElem)
  i = index_of_generator(x)
  
  return gen(base_algebra(x.parent), i)
end

function left_action_on_generator(g::MatrixGroupElem{T}, x::DrinfeldHeckeAlgebraElem{T, S}) where {T <: FieldElem, S <: RingElem}
  v = generator_to_vector(x)
  gv = matrix(g) * v
  
  return vector_to_algebra_element(x.parent, gv)
end



