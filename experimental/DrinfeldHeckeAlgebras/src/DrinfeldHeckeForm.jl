# TODO: Gibt es in Characteristic 0 ein Beispiel bei dem es nur die 0-DH form gibt?
# TODO: DOku schreiben
# TODO: Mehr Tests, z.b. über quadratic number fields, S5 über C
# TODO: Validierung ebenfalls über beide Strategien
# TODO: Code umstrukturieren, keine Validierung nach gnerischer Erstellung nötig

################################################################################
# Drinfeld-Hecke form
#
# Struct and methods for Drinfeld-Hecke form κ defined by 
#   κ = sum_g∈G (κ_g * g)
# where the κ_g are alternating bilinear forms fulfilling certain requirements
# (see Lemma 1.5 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002).
################################################################################

mutable struct DrinfeldHeckeForm{T <: FieldElem, S <: RingElem}
  group::MatrixGroup{T} # Matrix group G over a field K
  base_ring::Ring # Underlying K-algebra R with element type S of Drinfeld-Hecke form
  symmetric_algebra::MPolyRing{S} # Symmetric algebra R[V]
  group_algebra::GroupAlgebra # Group algebra R[V]G, used as base for R[V]#G
  forms::Dict{MatrixGroupElem{T}, BilinearForm{S}} # Alternating bilinear forms defining Drinfeld-Hecke form
  
  # Creates trivial (zero) Drinfeld-Hecke form from K-matrix-group
  function DrinfeldHeckeForm(G::MatrixGroup{T}) where {T <: FieldElem}
    return DrinfeldHeckeForm(G, base_ring(G))
  end

  # Creates trivial (zero) Drinfeld-Hecke form from K-matrix-group and K-algebra
  function DrinfeldHeckeForm(G::MatrixGroup{T}, R::Ring) where {T <: FieldElem}
    # Check if R is a K-algebra where K is the field over which G is defined
    K = base_ring(G)
    try
      # TODO
    catch
      throw(ArgumentError("The given ring must be an algebra over the base ring of the given group."))
    end
    
    S = elem_type(typeof(R))
    κ = new{T, S}()
        
    κ.group = G
    κ.base_ring = R
    κ.forms = Dict{MatrixGroupElem{T}, BilinearForm{S}}()
  
    RV, _ = polynomial_ring(R, ["x" * string(i) for i in 1:degree(G)])
    
    κ.symmetric_algebra = RV
    κ.group_algebra = RV[G]
    
    return κ
  end

  # Creates Drinfeld-Hecke form from non-empty forms input
  function DrinfeldHeckeForm(forms::Dict{MatrixGroupElem{T}, MatElem{S}}) where {T <: FieldElem, S <: RingElem}
    if length(forms) == 0
      throw(ArgumentError(
        "Forms must not be empty. To create zero form use 
        DrinfeldHeckeForm(G::MatrixGroup{T}) or DrinfeldHeckeForm(G::MatrixGroup{T}, R::Ring)"
      ))
    end

    g, κ_g = first(forms)
    G = g.parent
    R = base_ring(κ_g)
    
    κ = DrinfeldHeckeForm(G, R)
    set_forms(κ, forms)
  
    return κ
  end
end

################################################################################
# Drinfeld-Hecke form creation
################################################################################

const drinfeld_hecke_form = DrinfeldHeckeForm

function generic_drinfeld_hecke_form(G::MatrixGroup{T}) where {T <: RingElem}
  return generic_drinfeld_hecke_form(G, base_ring(G))
end

function generic_drinfeld_hecke_form(G::MatrixGroup{T}, R::Ring) where {T <: FieldElem}
  if R isa Field && characteristic(R) == 0
    forms = generate_generic_forms_locally(G, R)
  else
    forms = generate_generic_forms_globally(G, R)
  end

  if length(forms) == 0
    return DrinfeldHeckeForm(G, R)
  end
  
  return DrinfeldHeckeForm(forms)
end

################################################################################
# Drinfeld-Hecke form mutation
################################################################################

# Substitute parameters by values in any subring
function evaluate_parameters(κ::DrinfeldHeckeForm{T, S}, values::Vector) where {T <: FieldElem, S <: RingElem}
  if !(base_ring(κ) isa MPolyRing)
    throw(ArgumentError("Given form does not have any parameters."))
  end

  κ_evaluated = deepcopy_internal(κ)

  # If κ is zero, there is nothing to do
  if is_zero(κ) return κ_evaluated end
  
  R = base_ring(κ)
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
  forms = Dict{MatrixGroupElem{T}, MatElem{S}}()
  for (g, κ_g) in κ.forms
    forms[g] = map(λ, matrix(κ_g))
  end
  
  # Set them to κ_evaluated
  set_forms(κ_evaluated, forms)
  
  return κ_evaluated
end

# Checks if given forms define valid Drinfeld-Hecke form and if yes, sets them to κ
function set_forms(
  κ::DrinfeldHeckeForm{T, S}, 
  forms::Dict{MatrixGroupElem{T}, MatElem{S}}
) where {T <: FieldElem, S <: RingElem}
  # Check if forms input defines valid Drinfeld-Hecke form
  if !is_drinfeld_hecke_form(forms)
    throw(ArgumentError("The given alternating bilinear forms do not define a valid Drinfeld-Hecke form."))
  end
  
  # Initialize forms of κ
  κ.forms = Dict{MatrixGroupElem{T}, BilinearForm{S}}()
  
  # Set nonzero forms to κ
  for (g, m) in forms
    if !is_zero(m)
      κ.forms[g] = alternating_bilinear_form(m)
    end
  end
end

################################################################################
# String I/O
################################################################################

function show(io::IO, κ::DrinfeldHeckeForm)
  print(io, "Drinfeld-Hecke form given by ")

  if (length(κ.forms) == 0)
    print(io, "0")
    return
  end

  n = degree(group(κ))
  println(io, "alternating bilinear forms")
  for (_,(g, κ_g)) in enumerate(κ.forms)
    A = matrix(g)
    B = matrix(κ_g)
    
    for i in 1:n
      print(io, "   [")
      
      for j in 1:n
        mcl = max_column_length(A, j)
        print(io, repeat(" ", mcl - length(string(A[i,j]))))
        print(io, A[i,j])
        if j < n print(io, "   ") end
      end
      
      if i == n/2 || i == (n-1)/2
        print(io, "] => [")
      else
        print(io, "]    [")
      end
    
      for j in 1:n
        mcl = max_column_length(B, j)
        print(io, repeat(" ", mcl - length(string(B[i,j]))))
        print(io, B[i,j])
        if j < n print(io, "   ") end
      end
    
      print(io, "]")
      println()
    end
  
    println()
  end
end

function max_column_length(M::MatElem, j::Int)
  result = 0
  
  for i in 1:nrows(M)
    element_length = length(string(M[i,j]))
    
    if element_length > result
      result = element_length
    end
  end

  return result
end

################################################################################
# Generic functions
################################################################################

is_zero(κ::DrinfeldHeckeForm) = length(κ.forms) == 0
base_ring(κ::DrinfeldHeckeForm) = κ.base_ring
symmetric_algebra(κ::DrinfeldHeckeForm) = κ.symmetric_algebra
group(κ::DrinfeldHeckeForm) = κ.group
group_algebra(κ::DrinfeldHeckeForm) = κ.group_algebra
forms(κ::DrinfeldHeckeForm) = κ.forms
number_of_forms(κ::DrinfeldHeckeForm) = length(κ.forms)
nforms(κ::DrinfeldHeckeForm) = number_of_forms(κ)
parameters(κ::DrinfeldHeckeForm) = gens(base_ring(κ))

function deepcopy_internal(κ::DrinfeldHeckeForm)
  G = group(κ)
  R = base_ring(κ)
  
  κ_copy = DrinfeldHeckeForm(G, R)
  κ_copy.forms = κ.forms
  
  return κ_copy
end

function Base.getindex(κ::DrinfeldHeckeForm{T, S}, g::MatrixGroupElem{T}) where {T <: FieldElem, S <: RingElem}
  if haskey(κ.forms, g)
    return κ.forms[g]
  end

  d = degree(κ.group)
  m = zero_matrix(κ.base_ring, d, d)
  return alternating_bilinear_form(m)
end

function get_forms(κ::DrinfeldHeckeForm{T, S}) where {T <: FieldElem, S <: RingElem}
  forms = Dict{MatrixGroupElem{T}, MatElem{S}}()
  
  for (g, κ_g) in κ.forms
    forms[g] = matrix(κ_g)
  end

  return forms
end

################################################################################
# Application of Drinfeld-Hecke form
################################################################################

function (κ::DrinfeldHeckeForm{T, S})(v::Vector{S}, w::Vector{S}) where {T <: FieldElem, S <: RingElem}
  RG = κ.group_algebra
  RV = κ.symmetric_algebra
  
  res = RG()
  
  if v == w return res end

  for (g, κ_g) in κ.forms
    res = res + RV(κ_g(v,w)) * RG(g)
  end
  
  return res
end




