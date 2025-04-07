################################################################################
# Drinfeld Hecke Algebra
################################################################################

mutable struct DrinfeldHeckeAlgebra{T <: RingElem} <: NCRing
  drinfeld_hecke_form::DrinfeldHeckeForm{T}
  
  # Create from a Drinfeld-Hecke form
  function DrinfeldHeckeAlgebra(κ::DrinfeldHeckeForm{T}) where {T <: RingElem}
    return new{T}(κ)
  end

  # Create the zero Drinfeld-Hecke algebra from a group
  function DrinfeldHeckeAlgebra(G::MatrixGroup)
    κ = drinfeld_hecke_form(G)

    return DrinfeldHeckeAlgebra(κ)
  end
end

################################################################################
# Drinfeld Hecke Algebra Element
################################################################################

mutable struct DrinfeldHeckeAlgebraElem{T <: RingElem} <: NCRingElem
  parent::DrinfeldHeckeAlgebra{T}
  element::GroupAlgebraElem
  
  function DrinfeldHeckeAlgebraElem(A::DrinfeldHeckeAlgebra{T}, a::GroupAlgebraElem) where {T <: RingElem}
    if (a.parent != group_algebra(A))
        throw(ArgumentError("Element does not belong to the underlying group algebra")) 
    end

    return new{T}(A, a)
  end
end

################################################################################
# String I/O
################################################################################

function Base.show(io::IO, A::DrinfeldHeckeAlgebra)
  println(io, "Drinfeld-Hecke algebra")
  print(io, " of ")
  show(io, group(A))
  println(io)
  print(io, " defined by ")
  show(io, A.drinfeld_hecke_form)
end

function Base.show(io::IO, a::DrinfeldHeckeAlgebraElem)
  print(io, "Drinfeld-Hecke algebra element with coefficients ")
  show(io, a.element)
  # TODO nice repr. of element
end

################################################################################
# Generic functions
################################################################################

is_trivial(A::DrinfeldHeckeAlgebra) = is_zero(A.drinfeld_hecke_form)
dimension(A::DrinfeldHeckeAlgebra) = degree(group(A))
base_ring(A::DrinfeldHeckeAlgebra) = base_ring(A.drinfeld_hecke_form)
ring(A::DrinfeldHeckeAlgebra) = ring(A.drinfeld_hecke_form)
group(A::DrinfeldHeckeAlgebra) = group(A.drinfeld_hecke_form)
group_algebra(A::DrinfeldHeckeAlgebra) = group_algebra(A.drinfeld_hecke_form)

################################################################################
# Constructors
################################################################################

const drinfeld_hecke_algebra = DrinfeldHeckeAlgebra

function parametrized_drinfeld_hecke_algebra(G::MatrixGroup{T}) where {T <: RingElem}
  κ = parametrized_drinfeld_hecke_form(G)

  return drinfeld_hecke_algebra(κ)
end

(A::DrinfeldHeckeAlgebra)() = DrinfeldHeckeAlgebraElem(A, group_algebra(A)())
(A::DrinfeldHeckeAlgebra)(a::Integer) = DrinfeldHeckeAlgebraElem(A, group_algebra(A)(a))
(A::DrinfeldHeckeAlgebra)(a::RingElem) = DrinfeldHeckeAlgebraElem(A, group_algebra(A)(a))
(A::DrinfeldHeckeAlgebra)(a::MPolyRingElem) = DrinfeldHeckeAlgebraElem(A, group_algebra(A)(a))
(A::DrinfeldHeckeAlgebra)(g::MatrixGroupElem) = DrinfeldHeckeAlgebraElem(A, group_algebra(A)(g))
(A::DrinfeldHeckeAlgebra)(a::GroupAlgebraElem) = DrinfeldHeckeAlgebraElem(A, a)
function (A::DrinfeldHeckeAlgebra)(a::DrinfeldHeckeAlgebraElem)
  if a.parent == A 
    return a 
  else 
    throw(ArgumentError("Element does not belong to the given Drinfeld-Hecke algebra")) 
  end
end

################################################################################
# Data type and parent object methods
################################################################################

parent_type(::Type{DrinfeldHeckeAlgebraElem{T}}) where {T <: RingElem} = DrinfeldHeckeAlgebra{T}
elem_type(::Type{DrinfeldHeckeAlgebra{T}}) where {T <: RingElem} = DrinfeldHeckeAlgebraElem{T}
ring_type(::Type{DrinfeldHeckeAlgebra{T}}) where {T <: RingElem} = parent_type(T)
parent(a::DrinfeldHeckeAlgebraElem) = a.parent
is_domain_type(::Type{DrinfeldHeckeAlgebraElem}) = false
is_exact_type(::Type{DrinfeldHeckeAlgebraElem}) = true
Base.hash(a::DrinfeldHeckeAlgebraElem, h::UInt) = Base.hash(a.element, h)
deepcopy_internal(a::DrinfeldHeckeAlgebraElem, dict::IdDict) = DrinfeldHeckeAlgebraElem(a.parent, a.element)

################################################################################
# Basic manipulation of rings and elements
################################################################################

zero(A::DrinfeldHeckeAlgebra) = A()
one(A::DrinfeldHeckeAlgebra) = A(1)
iszero(a::DrinfeldHeckeAlgebraElem) = iszero(a.element)
isone(a::DrinfeldHeckeAlgebraElem) = isone(a.element)

################################################################################
# Canonicalisation
################################################################################

canonical_unit(a::DrinfeldHeckeAlgebraElem) = one(a.parent)

################################################################################
# Arithmetic Operations
################################################################################

-(a::DrinfeldHeckeAlgebraElem) = DrinfeldHeckeAlgebraElem(a.parent, -a.element)

+(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem) = DrinfeldHeckeAlgebraElem(a.parent, a.element + b.element)
-(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem) = DrinfeldHeckeAlgebraElem(a.parent, a.element - b.element)
*(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem) = multiply(a,b)

==(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem) = a.element == b.element
isequal(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem) = a == b

function ^(a::DrinfeldHeckeAlgebraElem, e::Int)
  if e == 0 return one(a.parent) end
  if e == 1 return a end
  
  return multiply(a,a^(e - 1))
end

################################################################################
# Unsafe operators
################################################################################

function zero!(a::DrinfeldHeckeAlgebraElem)
  a = zero(a.parent)
  return a
end

function mul!(c::DrinfeldHeckeAlgebraElem, a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem)
  c = a * b
  return c
end

function add!(c::DrinfeldHeckeAlgebraElem, a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem)
  c = a + b
  return c
end

function addeq!(a::DrinfeldHeckeAlgebraElem, b::DrinfeldHeckeAlgebraElem)
  a = a + b
  return a
end

################################################################################
# Random generation
################################################################################

rand(A::DrinfeldHeckeAlgebra) = DrinfeldHeckeAlgebraElem(A, rand(group_algebra(A)))

################################################################################
# Promotion rules
################################################################################

AbstractAlgebra.promote_rule(::Type{DrinfeldHeckeAlgebraElem{T}}, ::Type{DrinfeldHeckeAlgebraElem{T}}) where {T <: RingElem} = DrinfeldHeckeAlgebraElem{T}
function AbstractAlgebra.promote_rule(::Type{DrinfeldHeckeAlgebraElem{T}}, ::Type{U}) where {T <: RingElement, U <: RingElement}
 promote_rule(T, U) == T ? DrinfeldHeckeAlgebraElem{T} : Union{}
end
# TODO

################################################################################
# Multiplication in Drinfeld-Hecke Algebras
################################################################################

# Multiply group element g from the right
function *(a::DrinfeldHeckeAlgebraElem, g::MatrixGroupElem)
  A = a.parent
  
  return A(a.element * group_algebra(A)(g))
end

# Multiply with scalar from left and right
function *(c::T, a::DrinfeldHeckeAlgebraElem{T}) where {T <: RingElement}
  A = a.parent

  return A(ring(A)(c) * a.element)
end

*(a::DrinfeldHeckeAlgebraElem{T}, c::T) where {T <: RingElement} = c * a

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
  gx = generator_to_polynomial(x)^g
  
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
  κ = A.drinfeld_hecke_form
  xy = A(y.element * x.element) - κ(y,x)
  (c, m, g, tail) = split_element(xy)

  # m1 * m2 = mm1 * (x * y) * mm2 = mm1 * (c * m * g + tail) * mm2 = c * mm1 * m * g * mm2 + mm1 * tail * mm2
  return c * multiply(multiply_m1_with_m2(mm1, m) * g, mm2) + multiply(mm1, multiply_a_with_m(tail, mm2))
end

################################################################################
# Helper functions for recursion
################################################################################

# Returns quadruple (c, m, g, tail) where
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

################################################################################
# Helper functions for multiplication
################################################################################

(κ::DrinfeldHeckeForm)(x::DrinfeldHeckeAlgebraElem, y::DrinfeldHeckeAlgebraElem) = x.parent(κ(to_vector(x), to_vector(y)))

function to_vector(x::DrinfeldHeckeAlgebraElem)
  R = base_ring(x.parent)
  v = [R(0) for _ in 1:dimension(x.parent)]
  v[index_of_generator(x)] = R(1)
  
  return v
end

function index_of_generator(x::DrinfeldHeckeAlgebraElem)
  A = x.parent
  
  for i in 1:dimension(A)
    if A(gen(ring(A), i)) == x
      return i
    end
  end

  throw(ArgumentError("Element is not a generator."))
end

function generator_to_polynomial(x::DrinfeldHeckeAlgebraElem)
  i = index_of_generator(x)
  
  return gen(ring(x.parent), i)
end


