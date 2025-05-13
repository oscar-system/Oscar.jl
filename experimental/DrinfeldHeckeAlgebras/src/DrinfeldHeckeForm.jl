################################################################################
# Drinfeld-Hecke form
################################################################################

# Struct for Drinfeld-Hecke form κ defined by 
#   κ = sum_g∈G (κ_g * g)
# where the κ_g are alternating bilinear forms.
#
# Each Drinfeld-Hecke form defines a Drinfeld-Hecke algebra via
#   H_κ = R[V]#G / I_κ
# with
#   I_κ = < vw - wv - κ(v,w) | v,w ∈ V >
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
  forms = generate_generic_forms(G, R)

  if length(forms) == 0
    # There is only the trivial Drinfeld-Hecke form
    return DrinfeldHeckeForm(G, R)
  end
  
  return DrinfeldHeckeForm(forms)
end

################################################################################
# Drinfeld-Hecke form mutation
################################################################################

# Checks if given matrix defines a valid Drinfeld-Hecke form and if yes, sets it to κ and calculates according forms
# for the conjugates
function Base.setindex!(
  κ::DrinfeldHeckeForm{T, S}, 
  κ_g::MatElem{S}, 
  g::MatrixGroupElem{T}
) where {T <: FieldElem, S <: RingElem}
  if !is_valid_form(g, κ_g)
    throw(ArgumentError("Form does not define valid Drinfeld-Hecke form"))
  end
  
  for c in conjugacy_class(g)
    κ_c = calculate_form_for_conjugate(g, c, κ_g)
    
    if is_zero(κ_c)
      delete!(κ.forms, c)
    else 
      κ.forms[c] = alternating_bilinear_form(κ_c)
    end
  end
end

# Substitute parameters by values in any subring
function evaluate_parameters(κ::DrinfeldHeckeForm{T, S}, values::Vector) where {T <: FieldElem, S <: RingElem}
  if !is_generic(κ)
    throw(ArgumentError("Given form does not have any parameters."))
  end
  
  R = base_ring(κ)
  n = ngens(S)
  
  if length(values) != n
    throw(ArgumentError("Values input must contain exactly " * string(n) * " entries."))
  end
  
  # Check if values are in K
  safe_values = nothing
  try
    safe_values = map(v -> S(v), values)
  catch e
    throw(ArgumentError("The given values can not be cast into elements of the base ring."))
  end

  # Evaluation function
  φ = hom(S, S, safe_values)
  λ = f -> φ(f)

  # Apply homomorphism to forms
  forms = Dict{MatrixGroupElem{T}, MatElem{S}}()
  for (g, κ_g) in κ.forms
    forms[g] = map(λ, matrix(κ_g))
  end
  
  return DrinfeldHeckeForm(forms)
end

# Checks if given forms define valid Drinfeld-Hecke form and if yes, sets them to κ
function set_forms(
  κ::DrinfeldHeckeForm{T, S}, 
  forms::Dict{MatrixGroupElem{T}, MatElem{S}}
) where {T <: FieldElem, S <: RingElem}
  if !are_valid_forms(forms)
    throw(ArgumentError("Forms must define valid Drinfeld-Hecke form"))
  end
  
  # Set the forms to κ
  κ.forms = Dict{MatrixGroupElem{T}, BilinearForm{S}}()
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

  println(io, "alternating bilinear forms")
  for (i,(g, κ_g)) in enumerate(κ.forms)
    print(io, " ")
    print(io, g)
    print(io, " => ")
    println(io, κ_g.matrix)
  end
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
is_generic(κ::DrinfeldHeckeForm) = base_ring(κ) isa MPolyRing

function Base.getindex(κ::DrinfeldHeckeForm{T, S}, g::MatrixGroupElem{T}) where {T <: FieldElem, S <: RingElem}
  if haskey(κ.forms, g)
    return κ.forms[g]
  end

  d = degree(κ.group)
  m = zero_matrix(κ.base_ring, d, d)
  return alternating_bilinear_form(m)
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




