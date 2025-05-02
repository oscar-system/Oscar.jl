################################################################################
# Drinfeld-Hecke form
################################################################################

mutable struct DrinfeldHeckeForm{T <: RingElem}
  base_ring::Union{Ring, MPolyRing{T}}
  ring::MPolyRing
  group::MatrixGroup{T}
  group_algebra::GroupAlgebra
  forms::Dict{MatrixGroupElem, BilinearForm}
  values::Dict{Tuple{Vector,Vector}, GroupAlgebraElem}
  relation_handler::RelationHandler{T}

  # Creates Drinfeld-Hecke form from group, either zero or parametrized solution
  function DrinfeldHeckeForm(G::MatrixGroup{T}, parametrize::Bool = false) where {T <: RingElem}
    # Depending on the ring over which G is defined, GAP might not support iterating and computing inverses
    # In this case we throw an exception
    try
      elements(G)
      inv(G[1])
    catch
      throw(ArgumentError("Drinfeld-Hecke forms for groups over given ring are not supported."))
    end
  
    κ = new{T}()
    
    κ.relation_handler = RelationHandler(G)
    κ.group = G
    κ.forms = Dict{MatrixGroupElem, BilinearForm}()
    κ.values = Dict{Tuple{Vector,Vector}, GroupAlgebraElem}()
    
    d = degree(G)
    R = base_ring(G)
    
    if parametrize
      forms = solve_relations(κ.relation_handler)

      for (g, B) in forms
        if !is_zero(B)
          κ.forms[g] = alternating_bilinear_form(B)
        end
      end
    
      if length(forms) > 0
        _, B = first(forms)
        R = base_ring(B)
      end
    end
  
    κ.base_ring = R
  
    S, _ = polynomial_ring(R, ["x" * string(i) for i in 1:degree(G)])
    
    κ.ring = S
    κ.group_algebra = S[G]
    
    return κ
  end

  # Creates Drinfeld-Hecke form from non-empty forms input
  function DrinfeldHeckeForm(forms::Dict{MatrixGroupElem{T}, MatElem{T}}) where {T <: RingElem}
    if length(forms) == 0
      throw(ArgumentError("forms must not be empty. To create zero form use DrinfeldHeckeForm(G::MatrixGroup{T})")) 
    end
          
    g, _ = first(forms)
    G = g.parent
    
    κ = DrinfeldHeckeForm(G)
    
    set_forms(κ, forms)
  
    return κ
  end
end

################################################################################
# Drinfeld-Hecke form creation
################################################################################

const drinfeld_hecke_form = DrinfeldHeckeForm

function parametrized_drinfeld_hecke_form(G::MatrixGroup{T}) where {T <: RingElem}
  return DrinfeldHeckeForm(G, true)
end

################################################################################
# Parametrization helper functions
################################################################################

# Returns a parametrized solution to the linear equation system Mx = 0
function solve_relations(handler::RelationHandler{T}) where {T <: RingElem}
  M = relation_matrix(handler)
  G = group(handler)
  R = base_ring(M)
  n = ncols(M)
  K = zero_matrix(R, n, 0)
  
  try
    if is_zero(M)
      # We check this because over some rings nullspace is not implemented
      K = identity_matrix(R, n)
    else
      _, K = nullspace(M)
    end
  catch e
    throw(ArgumentError("Can not solve Drinfeld-Hecke relations over given ring."))
  end

  n = ncols(K)
  m = nrows(K)
  forms = Dict{MatrixGroupElem, MatElem}()
  
  # If n is zero, we don't need parameters since zero is the only solution
  if n == 0
    return forms
  end
        
  # For creating a parametrized solution, we work over a polynomial ring S = R[t1,...tn]
  parameters = n == 1 ? ["t"] : ["t" * string(i) for i in 1:n]
  S, _ = polynomial_ring(R, parameters)

  # Use kernel basis to create a parametrized solution for a Drinfeld-Hecke form
  x = fill(S(), m)
  for j in 1:n
    x = x + [S[j] * K[i, j] for i in 1:m]
  end

  # Initialize forms
  d = degree(G)
  for g in G
    forms[g] = zero_matrix(S,d,d)
  end
    
  # Use handler map to convert result to matrices
  for ((g,i,j),k) in map(handler)
    forms[g][i,j] = x[k]
    forms[g][j,i] = -x[k]
  end
  
  return forms 
end

################################################################################
# Drinfeld-Hecke form mutation
################################################################################

function Base.setindex!(κ::DrinfeldHeckeForm{T}, m::MatElem{T}, g::MatrixGroupElem{T}) where {T <: RingElem}
  forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
  
  for (h, κ_h) in κ.forms
    forms[h] = matrix(κ_h)
  end

  forms[g] = m
  
  set_forms(κ, forms)
end

function evaluate_parameters(κ::DrinfeldHeckeForm{T}, values::Vector) where {T <: RingElem}
  if !is_parametrized(κ)
    throw(ArgumentError("Given form does not have any parameter."))
  end
  
  S = base_ring(κ)
  n = ngens(S)
  
  if length(values) != n
    throw(ArgumentError("Values input must contain exactly " * string(n) * " entries."))
  end

  R = base_ring(S)
  
  # Check if values are in K
  safe_values = nothing
  try
    safe_values = map(v -> R(v), values)
  catch e
    throw(ArgumentError("The given values can not be cast into elements of the base ring."))
  end

  # Evaluation function
  φ = hom(S, R, safe_values)
  λ = f -> φ(f)

  # Move forms to new ring S_new
  forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
  for (g, κ_g) in κ.forms
    forms[g] = map(λ, matrix(κ_g))
  end
  
  return DrinfeldHeckeForm(forms)
end

function set_forms(κ::DrinfeldHeckeForm{T}, forms::Dict{MatrixGroupElem{T}, MatElem{T}}) where {T <: RingElem}
  if is_parametrized(κ)
    throw(ArgumentError("Can not set forms to parametrized Drinfeld-Hecke form, please use evaluate_parameters"))
  end
  
  G = κ.group
  R = κ.base_ring
  map = κ.relation_handler.map
  M = κ.relation_handler.relation_matrix
  
  # We need to check if forms satisfy the needed relations to define a Drinfeld-Hecke form
  # For this we translate the forms into a vector x using the handler map
  x = fill(R(), length(map))
  for ((g,i,j),k) in map
    if haskey(forms, g)
      x[k] = forms[g][i,j]
    end
  end
  
  # Check if x satisfies the relations
  if !is_zero(M*x)
    throw(ArgumentError("The given forms do not define a valid Drinfeld-Hecke form"))
  end
  
  # Set the forms to κ
  κ.forms = Dict{MatrixGroupElem{T}, BilinearForm{T}}()
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
ring(κ::DrinfeldHeckeForm) = κ.ring
group(κ::DrinfeldHeckeForm) = κ.group
group_algebra(κ::DrinfeldHeckeForm) = κ.group_algebra
forms(κ::DrinfeldHeckeForm) = κ.forms
number_of_forms(κ::DrinfeldHeckeForm) = length(κ.forms)
nforms(κ::DrinfeldHeckeForm) = number_of_forms(κ)

function is_parametrized(κ::DrinfeldHeckeForm{T}) where {T <: RingElem}
  R = base_ring(κ)
  G = group(κ)
  
  return R isa MPolyRing{T} && G isa MatrixGroup{T}
end

function parameters(κ::DrinfeldHeckeForm)
  if !is_parametrized(κ)
    throw(ArgumentError("Given form does not have any parameter."))
  end
  
  S = base_ring(κ)
  return gens(S)
end

function Base.getindex(κ::DrinfeldHeckeForm{T}, g::MatrixGroupElem{T}) where {T <: RingElem}
  if haskey(κ.forms, g) 
    return κ.forms[g] 
  else 
    d = degree(κ.group)
    m = zero_matrix(κ.base_ring, d, d)
    return alternating_bilinear_form(m)
  end
end

################################################################################
# Application of Drinfeld-Hecke form
################################################################################

function (κ::DrinfeldHeckeForm{T})(v::Vector{T}, w::Vector{T}) where T <: RingElem
  RG = κ.group_algebra
  
  if v == w return RG() end
  
  if haskey(κ.values, (v,w))
    return κ.values[(v,w)]
  end

  if haskey(κ.values, (w,v))
    κ.values[(v,w)] = -1 * κ.values[(w,v)]
    return κ.values[(v,w)]
  end
  
  res = RG()

  for (g, κ_g) in κ.forms
    res = res + base_ring(RG)(κ_g(v,w)) * RG(g)
  end

  κ.values[(v,w)] = res
  
  return res
end




