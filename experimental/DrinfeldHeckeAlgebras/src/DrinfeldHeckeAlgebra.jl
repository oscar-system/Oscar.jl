################################################################################
# Struct and methods to generate and handle concrete or generic Drinfeld-Hecke algebras
#
# Cassandra Koenen, 2025
################################################################################

################################################################################
# Struct and different constructors for Drinfeld-Hecke algebras
################################################################################

@doc raw"""
    DrinfeldHeckeAlgebra{T <: FieldElem, S <: RingElem} <: NCRing

Type for representing a Drinfeld-Hecke algebra for a group $G$. The parameter ```T``` is the element type over which
$G$ is defined, ```S``` is the ring element type over which the associated Drinfeld-Hecke form is defined.
"""
mutable struct DrinfeldHeckeAlgebra{T <: FieldElem, S <: RingElem} <: NCRing
  form::DrinfeldHeckeForm{T}
  
  # Create from forms input
  function DrinfeldHeckeAlgebra(forms::Dict)
    κ = DrinfeldHeckeForm(forms)

    return DrinfeldHeckeAlgebra(κ)
  end

  # Create the zero Drinfeld-Hecke algebra from a group and a ring
  function DrinfeldHeckeAlgebra(G::MatrixGroup{T}, R::Ring=base_ring(G)) where {T <: FieldElem}
    κ = DrinfeldHeckeForm(G, R)

    return DrinfeldHeckeAlgebra(κ)
  end

  # Create from a Drinfeld-Hecke form
  function DrinfeldHeckeAlgebra(κ::DrinfeldHeckeForm{T, S}) where {T <: FieldElem, S <: RingElem}
    return new{T, S}(κ)
  end
end

@doc raw"""
    drinfeld_hecke_algebra(G::MatrixGroup{T}, R::Ring=base_ring(G)) where {T <: FieldElem}

Create the trivial Drinfeld-Hecke form for the matrix group ```G``` over the ring ```R```, i.e. the skew group ring $R[V]\#G$ 
where $V$ is the vector space on which ```G``` acts

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> R, _ = polynomial_ring(QQ, ["x","y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> drinfeld_hecke_algebra(G, R)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   x, y
given by 0
```
"""
drinfeld_hecke_algebra(G::MatrixGroup{T}, R::Ring=base_ring(G)) where {T <: FieldElem} = DrinfeldHeckeAlgebra(G, R)

@doc raw"""
    drinfeld_hecke_algebra(forms::Dict)

Create the Drinfeld-Hecke algebra given by the dictionary ```forms``` which is expected to have keys of type 
```MatrixGroupElem{T}``` and values of type ```MatElem{S}``` where ```T <: FieldElem, S <: RingElem```. If the types
are not correct, an exception is thrown.

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> R, (x,y) = polynomial_ring(QQ, ["x","y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> κ_1 = matrix(R, [0 x; -x 0])
[ 0   x]
[-x   0]

julia> κ_g = matrix(R, [0 y; -y 0])
[ 0   y]
[-y   0]

julia> forms = Dict(one(G) => κ_1, G[1] => κ_g)
Dict{MatrixGroupElem{QQFieldElem, QQMatrix}, AbstractAlgebra.Generic.MatSpaceElem{QQMPolyRingElem}} with 2 entries:
  [1 0; 0 1]   => [0 x; -x 0]
  [-1 0; 0 -1] => [0 y; -y 0]

julia> drinfeld_hecke_algebra(forms)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   x, y
given by alternating bilinear forms
   [1   0] => [ 0   x]
   [0   1]    [-x   0]

   [-1    0] => [ 0   y]
   [ 0   -1]    [-y   0]
```
"""
drinfeld_hecke_algebra(forms::Dict) = DrinfeldHeckeAlgebra(forms)

################################################################################
# Drinfeld Hecke Algebra Element
################################################################################

@doc raw"""
    DrinfeldHeckeAlgebraElem{T <: FieldElem, S <: RingElem} <: NCRingElem

Type for representing a Drinfeld-Hecke algebra element.
"""
mutable struct DrinfeldHeckeAlgebraElem{T <: FieldElem, S <: RingElem} <: NCRingElem
  parent::DrinfeldHeckeAlgebra{T, S}
  element::GroupAlgebraElem
  
  function DrinfeldHeckeAlgebraElem(A::DrinfeldHeckeAlgebra{T, S}, a::GroupAlgebraElem) where {T <: FieldElem, S <: RingElem}
    if (a.parent != group_algebra(A))
        throw(ArgumentError("Element does not belong to the underlying group algebra")) 
    end

    return new{T, S}(A, a)
  end
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
# Generic Drinfeld-Hecke algebra generation
################################################################################

@doc raw"""
    generic_drinfeld_hecke_algebra(G::MatrixGroup{T}) where {T <: FieldElem}

Create a generic (parametrized) Drinfeld-Hecke algebra for the group ```G```.

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> generic_drinfeld_hecke_algebra(G)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   t1, t2
given by alternating bilinear forms
   [1   0] => [  0   t1]
   [0   1]    [-t1    0]

   [-1    0] => [  0   t2]
   [ 0   -1]    [-t2    0]
```
"""
function generic_drinfeld_hecke_algebra(G::MatrixGroup{T}, R::Ring=base_ring(G)) where {T <: FieldElem}
  κ = generic_drinfeld_hecke_form(G,R)

  return DrinfeldHeckeAlgebra(κ)
end

################################################################################
# Evaluate parameters of generic Drinfeld-Hecke algebra
################################################################################

@doc raw"""
    evaluate_parameters(A::DrinfeldHeckeAlgebra, values::Vector)

Evaluate the parameters of the Drinfeld-Hecke algebra ```A``` by the values given in ```values```.

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> A = generic_drinfeld_hecke_algebra(G)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   t1, t2
given by alternating bilinear forms
   [1   0] => [  0   t1]
   [0   1]    [-t1    0]

   [-1    0] => [  0   t2]
   [ 0   -1]    [-t2    0]

julia> evaluate_parameters(A, [1,2])
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   t1, t2
given by alternating bilinear forms
   [1   0] => [ 0   1]
   [0   1]    [-1   0]

   [-1    0] => [ 0   2]
   [ 0   -1]    [-2   0]
```
"""
function evaluate_parameters(A::DrinfeldHeckeAlgebra, values::Vector)
  G = group(A)
  R = base_ring(A)
  
  if !(R isa MPolyRing)
    throw(ArgumentError("The given form does not have any parameters."))
  end

  # If A is zero, there is nothing to do
  if is_zero(form(A)) return drinfeld_hecke_algebra(G,R) end
  
  n = ngens(R)
  
  if length(values) != n
    throw(ArgumentError("Values input must contain exactly " * string(n) * " entries."))
  end
  
  # Check if values are in R
  safe_values = nothing
  try
    safe_values = map(v -> R(v), values)
  catch e
    throw(ArgumentError("The given values can not be cast into elements of the base ring."))
  end

  # Evaluation function
  φ = hom(R, R, safe_values)
  λ = f -> φ(f)

  # Apply homomorphism to forms
  forms = Dict()
  for (g, κ_g) in alternating_bilinear_forms(form(A))
    forms[g] = map(λ, matrix(κ_g))
  end
  
  # Create new DH algebra from forms
  return DrinfeldHeckeAlgebra(forms)
end

################################################################################
# String I/O
################################################################################

function Base.show(io::IO, A::DrinfeldHeckeAlgebra)
  println(io, "Drinfeld-Hecke algebra")
  print(io, "   for ")
  show(io, group(A))
  println(io)
  println(io, "with generators")
  print(io, "   " * join(gens(base_algebra(A)), ", "))
  println(io, ", " * join(["g" * string(i) for i in 1:ngens(group(A))], ", "))
  println(io)
  print(io, "defined by ")
  show(io, form(A))
end

function Base.show(io::IO, a::DrinfeldHeckeAlgebraElem)
  if is_zero(a)
    print(io, "0")
    return
  end

  if is_one(a)
    print(io, "1")
    return
  end

  group_elements = parent(a.element).base_to_group
  non_zero_coefficients = Dict()
  
  for i in 1:length(coefficients(a.element))
    c = coefficients(a.element)[i]
    
    if !is_zero(c)
      g = group_elements[i]
      non_zero_coefficients[g] = c
    end
  end

  for (i,(g,c)) in enumerate(non_zero_coefficients)
    if !is_one(c)
      if length(terms(c)) > 1
        print(io, "(")
      end
    
      print(io, c)
      
      if length(terms(c)) > 1
        print(io, ")")
      end
    end
  
    if !is_one(c) && !is_one(g)
      print(io, " * ")
    end
  
    if !is_one(g)
      print(io, g)
    end
    
    if i < length(non_zero_coefficients)
      print(io, " + ")
    end
  end
end

################################################################################
# Basic functionality
################################################################################

@doc raw"""
    base_ring(A::DrinfeldHeckeAlgebra)

Return the base ring over which ```A``` is defined.

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> R, (x,y) = polynomial_ring(QQ, ["x","y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> A = drinfeld_hecke_algebra(G, R)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   x, y
given by 0

julia> base_ring(A)
Multivariate polynomial ring in 2 variables x, y
  over rational field
```
"""
base_ring(A::DrinfeldHeckeAlgebra) = base_ring(form(A))

@doc raw"""
    base_algebra(A::DrinfeldHeckeAlgebra)

Return the polynomial ring over the base ring of ```A``` with indeterminants a standard 
basis of $V$, the vector space on which the group of ```A``` acts on.

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> A = drinfeld_hecke_algebra(G)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Rational field
given by 0

julia> base_algebra(A)
Multivariate polynomial ring in 2 variables x1, x2
  over rational field
```
"""
base_algebra(A::DrinfeldHeckeAlgebra) = base_algebra(form(A))

@doc raw"""
    group(A::DrinfeldHeckeAlgebra)

Return the group for which ```A``` is defined.

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> A = drinfeld_hecke_algebra(G)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Rational field
given by 0

julia> group(A)
Matrix group of degree 2
  over rational field
```
"""
group(A::DrinfeldHeckeAlgebra) = group(form(A))

group_algebra(A::DrinfeldHeckeAlgebra) = group_algebra(form(A))
form(A::DrinfeldHeckeAlgebra) = A.form

@doc raw"""
    generators(A::DrinfeldHeckeAlgebra)

Return the generators of ```A```.

Alias: ```gens```

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> A = drinfeld_hecke_algebra(G)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Rational field
given by 0

julia> generators(A)
3-element Vector{Oscar.DrinfeldHeckeAlgebraElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 [-1 0; 0 -1]
```
"""
function generators(A::DrinfeldHeckeAlgebra)
  indeterminants = map(x -> A(x), gens(base_algebra(A)))
  group_generators = map(x -> A(x), gens(group(A)))
  
  return vcat(indeterminants, group_generators)
end

gens(A::DrinfeldHeckeAlgebra) = generators(A)

function getindex(A::DrinfeldHeckeAlgebra, i::Int)
  return gens(A)[i]
end

@doc raw"""
    parameters(A::DrinfeldHeckeAlgebra)

Return the parameters of ```A```.

Alias: ```params```

# Examples
```jldoctest
julia> G = matrix_group(matrix(QQ, [-1 0;0 -1]))
Matrix group of degree 2
  over rational field

julia> A = generic_drinfeld_hecke_algebra(G)
Drinfeld-Hecke algebra
   for Matrix group of degree 2 over QQ
with generators
   x1, x2, g1

defined by Drinfeld-Hecke form over base ring
   Multivariate polynomial ring in 2 variables over QQ
with parameters 
   t1, t2
given by alternating bilinear forms
   [1   0] => [  0   t1]
   [0   1]    [-t1    0]

   [-1    0] => [  0   t2]
   [ 0   -1]    [-t2    0]

julia> parameters(A)
2-element Vector{QQMPolyRingElem}:
 t1
 t2
```
"""
parameters(A::DrinfeldHeckeAlgebra) = parameters(form(A))
params(A::DrinfeldHeckeAlgebra) = parameters(A)

################################################################################
# Data type and parent object methods
################################################################################

parent_type(::Type{DrinfeldHeckeAlgebraElem{T, S}}) where {T <: FieldElem, S <: RingElem} = DrinfeldHeckeAlgebra{T, S}
elem_type(::Type{DrinfeldHeckeAlgebra{T, S}}) where {T <: FieldElem, S <: RingElem} = DrinfeldHeckeAlgebraElem{T, S}
ring_type(::Type{DrinfeldHeckeAlgebra{T, S}}) where {T <: FieldElem, S <: RingElem} = parent_type(T)
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

AbstractAlgebra.promote_rule(::Type{DrinfeldHeckeAlgebraElem{T, S}}, ::Type{DrinfeldHeckeAlgebraElem{T, S}}) where {T <: FieldElem, S <: RingElem} = DrinfeldHeckeAlgebraElem{T, S}
function AbstractAlgebra.promote_rule(::Type{DrinfeldHeckeAlgebraElem{T, S}}, ::Type{U}) where {T <: FieldElem, S <: RingElem, U <: RingElement}
 promote_rule(T, U) == T ? DrinfeldHeckeAlgebraElem{T, S} : Union{}
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
function *(c::S, a::DrinfeldHeckeAlgebraElem{T, S}) where {T <: FieldElem, S <: RingElem}
  A = a.parent

  return A(base_algebra(A)(c) * a.element)
end

*(a::DrinfeldHeckeAlgebraElem{T, S}, c::S) where {T <: FieldElem, S <: RingElem} = c * a

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
  gx = x^g
  
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

function ^(x::DrinfeldHeckeAlgebraElem{T, S}, g::MatrixGroupElem{T}) where {T <: FieldElem, S <: RingElem}
  v = generator_to_vector(x)
  gv = matrix(g) * v
  
  return vector_to_algebra_element(x.parent, gv)
end



