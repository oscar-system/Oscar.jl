################################################################################
# Drinfeld-Hecke form
################################################################################

mutable struct DrinfeldHeckeForm{T <: RingElem}
  base_ring::Ring
  ring::MPolyRing{T}
  group::MatrixGroup{T}
  group_algebra::GroupAlgebra
  forms::Dict{MatrixGroupElem{T}, BilinearForm{T}}
  values::Dict{Tuple{Vector{T},Vector{T}}, GroupAlgebraElem}
  relation_handler::RelationHandler{T}

  # Creates empty (zero) Drinfeld-Hecke form
  function DrinfeldHeckeForm(G::MatrixGroup{T}) where {T <: RingElem}
    R = base_ring(G)
    d = degree(G)
    S, _ = polynomial_ring(R, ["x" * string(i) for i in 1:d])

    κ = new{T}()

    κ.base_ring = R
    κ.ring = S
    κ.group = G
    κ.group_algebra = S[G]
    κ.forms = Dict{MatrixGroupElem{T}, BilinearForm{T}}()
    κ.values = Dict{Tuple{Vector{T},Vector{T}}, GroupAlgebraElem}()
    κ.relation_handler = RelationHandler(G)
    
    return κ
  end

  # Creates concrete Drinfeld-Hecke form from non-empty forms input
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
  handler = RelationHandler(G)
  
  G = handler.group
  R = base_ring(G)
  d = degree(G)
  
  # Try to solve relations
  n, kernel = nullspace(handler.relation_matrix) # TODO what if R is not a field?
  
  # There is no solution except zero
  if n == 0 
    return DrinfeldHeckeForm(G) 
  end

  # For creating a parametrized solution, we work over a polynomial ring S = R[t1,...tn]
  S, _ = polynomial_ring(R, ["t" * string(i) for i in 1:n])
  
  # Use kernel basis to create a parametrized solution
  m = nrows(kernel)
  sol = fill(S(), m)
  for j in 1:n
    sol = sol + [S[j] * kernel[i, j] for i in 1:m]
  end

  # Cast group to new base ring S
  G_S = change_base_ring(S, G)
  T_S = elem_type(typeof(S))
  
  # Initialize forms
  forms = Dict{MatrixGroupElem{T_S}, MatElem{T_S}}()
  for g_S in G_S
    forms[g_S] = zero_matrix(S,d,d)
  end
    
  # Use handler map to convert result to matrices
  for ((g,i,j),k) in handler.map
    g_S = G_S((x -> S(x)).(matrix(g))) # cast group element g in G to an element of G_S
    forms[g_S][i,j] = sol[k]
    forms[g_S][j,i] = -sol[k]
  end
  
  return DrinfeldHeckeForm(forms)
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

function Base.getindex(κ::DrinfeldHeckeForm{T}, g::MatrixGroupElem{T}) where {T <: RingElem}
  if haskey(κ.forms, g) 
    return κ.forms[g] 
  else 
    d = degree(κ.group)
    m = zero_matrix(κ.base_ring, d, d)
    return alternating_bilinear_form(m)
  end
end

function Base.setindex!(κ::DrinfeldHeckeForm{T}, m::MatElem{T}, g::MatrixGroupElem{T}) where {T <: RingElem}
  forms = Dict{MatrixGroupElem{T}, MatElem{T}}()
  
  for (h, κ_h) in κ.forms
    forms[h] = matrix(κ_h)
  end

  forms[g] = m
  
  set_forms(κ, forms)
end

function set_forms(κ::DrinfeldHeckeForm{T}, forms::Dict{MatrixGroupElem{T}, MatElem{T}}) where {T <: RingElem}
  G = κ.group
  R = κ.base_ring
  map = κ.relation_handler.map
  M = κ.relation_handler.relation_matrix
  
  # We need to check if forms satisfies the needed relations to define a Drinfeld-Hecke form
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




