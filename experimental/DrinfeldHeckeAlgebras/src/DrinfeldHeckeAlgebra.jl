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

(A::DrinfeldHeckeAlgebra)() = A(group_algebra(A)())
(A::DrinfeldHeckeAlgebra)(a::Integer) = A(group_algebra(A)(base_algebra(A)(a)))
(A::DrinfeldHeckeAlgebra{T, S})(a::T) where {T <: FieldElem, S <: RingElem} = A(group_algebra(A)(base_algebra(A)(a)))
(A::DrinfeldHeckeAlgebra{T, S})(a::S) where {T <: FieldElem, S <: RingElem} = A(group_algebra(A)(base_algebra(A)(a)))
(A::DrinfeldHeckeAlgebra{T, S})(a::MPolyRingElem{S}) where {T <: FieldElem, S <: RingElem} = A(group_algebra(A)(a))
(A::DrinfeldHeckeAlgebra{T, S})(g::MatrixGroupElem{T}) where {T <: FieldElem, S <: RingElem} = A(group_algebra(A)(g))
(A::DrinfeldHeckeAlgebra)(a::GroupAlgebraElem) = DrinfeldHeckeAlgebraElem(A, a)

function (A::DrinfeldHeckeAlgebra)(a::DrinfeldHeckeAlgebraElem)
  if a.parent == A 
    return a 
  else 
    throw(ArgumentError("Element does not belong to the given Drinfeld-Hecke algebra")) 
  end
end

function (A::DrinfeldHeckeAlgebra)(a)
  A(group_algebra(A)(a))
  try
    return A(group_algebra(A)(a))
  catch 
    throw(ArgumentError("Element cannot be cast into the given Drinfeld-Hecke algebra")) 
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

  n = degree(group(parent(a)))
  for (i,(g,c)) in enumerate(non_zero_coefficients)
    if is_one(g)
      println(io, c)
    else
      # Print row by row
      for i in 1:n
        # Print c only if c is not one
        if !is_one(c)
          # Calculate how much space it takes up to print "c * "
          needed_space_to_print_c = 3 + length(string(c))
          if length(terms(c)) > 1
            needed_space_to_print_c += 2
          end
          
          if (i == n/2 || i == (n-1)/2)
            # Print c
            if length(terms(c)) > 1
              print(io, "(")
            end
          
            print(io, c)
            
            if length(terms(c)) > 1
              print(io, ") * ")
            else
              print(io, " * ")
            end
          else
            # Otherwise put whitespace
            print(io, repeat(" ", needed_space_to_print_c))
          end
        end
        
        if !is_one(g)
          # Start printing current row of g
          print(io, "[")
          A = matrix(g)
        
          for j in 1:n
            mcl = max_column_length(A, j)
            print(io, repeat(" ", mcl - length(string(A[i,j]))))
            print(io, A[i,j])
            if j < n print(io, "   ") end
          end
        
          println(io, "]")
        else
          println(io)
        end
      end
    end

    # print "+"
    if i < length(non_zero_coefficients)
      println(io, " + ")
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

function getindex(A::DrinfeldHeckeAlgebra, i::Int)
  generators = vcat(gens(base_algebra(A)), gens(group(A)))
  
  return A(generators[i])
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

@doc raw"""
    is_trivial(A::DrinfeldHeckeAlgebra)

Return ```true``` if ```A``` is the trivial Drinfeld-Hecke algebra, ```false``` otherwise.

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
is_trivial(A::DrinfeldHeckeAlgebra) = is_zero(form(A))

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
# Type coercion
################################################################################

*(a::Union{RingElem, MatrixGroupElem}, b::DrinfeldHeckeAlgebraElem) = b.parent(a) * b.parent(b)
*(a::DrinfeldHeckeAlgebraElem, b::Union{RingElem, MatrixGroupElem}) = a.parent(a) * a.parent(b)
+(a::Union{RingElem, MatrixGroupElem}, b::DrinfeldHeckeAlgebraElem) = b.parent(a) + b.parent(b)
+(a::DrinfeldHeckeAlgebraElem, b::Union{RingElem, MatrixGroupElem}) = a.parent(a) + a.parent(b)
-(a::Union{RingElem, MatrixGroupElem}, b::DrinfeldHeckeAlgebraElem) = b.parent(a) - b.parent(b)
-(a::DrinfeldHeckeAlgebraElem, b::Union{RingElem, MatrixGroupElem}) = a.parent(a) - a.parent(b)

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




