#######################################
# Struct and methods to generate and handle concrete or generic Drinfeld-Hecke forms
#
# Cassandra Koenen, 2025
#######################################

#######################################
# Struct and different constructors for Drinfeld-Hecke forms
#######################################

mutable struct DrinfeldHeckeForm{T <: FieldElem, S <: RingElem}
  group::MatrixGroup{T} # Matrix group G over a field K
  base_ring::Ring # Underlying K-algebra R with element type S of Drinfeld-Hecke form
  base_algebra::MPolyRing{S} # Polynomial ring R[V]
  group_algebra::GroupAlgebra # Group algebra R[V]G, used as base for R[V]#G
  forms::Dict{MatrixGroupElem{T}, AlternatingBilinearForm{S}} # Alternating bilinear forms defining Drinfeld-Hecke form

  # Creates trivial (zero) Drinfeld-Hecke form from K-matrix-group and K-algebra
  function DrinfeldHeckeForm(G::MatrixGroup{T}, R::Ring=base_ring(G)) where {T <: FieldElem}
    # Check if R contains K where K is the field over which G is defined
    K = base_ring(G)
    try
      R(one(K))
    catch
      throw(ArgumentError("The given ring must be an algebra over the base ring of the given group."))
    end
    
    S = elem_type(typeof(R))
    kappa = new{T, S}()
        
    kappa.group = G
    kappa.base_ring = R
    kappa.forms = Dict{MatrixGroupElem{T}, AlternatingBilinearForm{S}}()
  
    RV, _ = polynomial_ring(R, ["x" * string(i) for i in 1:degree(G)])
    
    kappa.base_algebra = RV
    kappa.group_algebra = RV[G]
    
    return kappa
  end

  # Creates Drinfeld-Hecke form from non-empty forms input
  function DrinfeldHeckeForm(forms::Dict)
    if length(forms) == 0
      throw(ArgumentError(
        "Forms must not be empty. To create zero form use 
        DrinfeldHeckeForm(G::MatrixGroup{T}) or DrinfeldHeckeForm(G::MatrixGroup{T}, R::Ring)"
      ))
    end
  
    g, kappa_g = first(forms)
    G = parent(g)
    R = base_ring(kappa_g)
    
    T = elem_type(typeof(base_ring(G)))
    S = elem_type(typeof(R))
    
    forms_safe = Dict{MatrixGroupElem{T}, MatElem{S}}()
    
    for (g,m) in forms
      if !(g isa MatrixGroupElem{T}) || !(m isa MatElem{S})
        throw(ArgumentError(
          "The forms dictionary must be of the structure 
          MatrixGroupElem{T} => MatElem{S} where {T <: FieldElem, S <: RingElem}"
        ))
      end
    
      forms_safe[g] = m
    end
  
    # Check if forms input defines valid Drinfeld-Hecke form
    if !is_drinfeld_hecke_form(forms_safe)
      throw(ArgumentError("The given forms do not define a valid Drinfeld-Hecke form."))
    end
  
    kappa = DrinfeldHeckeForm(G, R)
    
    for (g, m) in forms_safe
      if !is_zero(m)
        kappa.forms[g] = alternating_bilinear_form(m)
      end
    end
  
    return kappa
  end
end

#######################################
# Generic Drinfeld-Hecke form creation
#######################################

function generic_drinfeld_hecke_form(G::MatrixGroup{T}, R::Ring=base_ring(G)) where {T <: FieldElem}
  if R isa Field && characteristic(R) == 0
    forms = generate_generic_forms_locally(G, R)
  else
    forms = generate_generic_forms_globally(G, R)
  end

  # Filter nonzero forms
  for (g, kappa_g) in forms
    if is_zero(kappa_g)
      delete!(forms, g)
    end
  end

  # If the forms are empty, return zero algebra
  if length(forms) == 0
    return DrinfeldHeckeForm(G, R)
  end

  # Otherwise extract ring data
  g, kappa_g = first(forms)
  S = base_ring(kappa_g)
  
  # Create form and set forms
  kappa = DrinfeldHeckeForm(forms)
  
  return kappa
end

#######################################
# String I/O
#######################################

function show(io::IO, kappa::DrinfeldHeckeForm)
  println(io, "Drinfeld-Hecke form over base ring")
  println(io, "   " * string(base_ring(kappa)))

  if number_of_parameters(kappa) > 0
    println(io, "with parameters ")
    println(io, "   " * join(parameters(kappa), ", "))
  end

  if (length(kappa.forms) == 0)
    print(io, "given by 0")
    return
  end

  println(io, "given by alternating bilinear forms")
  n = degree(group(kappa))
  for (k,(g, kappa_g)) in enumerate(kappa.forms)
    A = matrix(g)
    B = matrix(kappa_g)
    
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
    
      if i == n && k == length(kappa.forms)
        print(io, "]")
      else
        println(io, "]")
      end
    end
  
    if k < length(kappa.forms)
      println(io)
    end
  end
end

#######################################
# Generic functions
#######################################

is_zero(kappa::DrinfeldHeckeForm) = length(kappa.forms) == 0
base_field(kappa::DrinfeldHeckeForm) = base_ring(kappa.group)
base_ring(kappa::DrinfeldHeckeForm) = kappa.base_ring
base_algebra(kappa::DrinfeldHeckeForm) = kappa.base_algebra
group(kappa::DrinfeldHeckeForm) = kappa.group
group_algebra(kappa::DrinfeldHeckeForm) = kappa.group_algebra
alternating_bilinear_forms(kappa::DrinfeldHeckeForm) = kappa.forms
number_of_forms(kappa::DrinfeldHeckeForm) = length(kappa.forms)
nforms(kappa::DrinfeldHeckeForm) = number_of_forms(kappa)
parameters(kappa::DrinfeldHeckeForm) = if base_ring(kappa) isa MPolyRing return gens(base_ring(kappa)) else return [] end
number_of_parameters(kappa::DrinfeldHeckeForm) = length(parameters(kappa))
nparams(kappa::DrinfeldHeckeForm) = number_of_parameters(kappa)

#######################################
# Application of Drinfeld-Hecke form
#######################################

function (kappa::DrinfeldHeckeForm{T, S})(v::Vector{S}, w::Vector{S}) where {T <: FieldElem, S <: RingElem}
  RG = group_algebra(kappa)
  RV = base_algebra(kappa)
  
  res = RG()
  
  if v == w return res end

  for (g, kappa_g) in alternating_bilinear_forms(kappa)
    res = res + RV(kappa_g(v,w)) * RG(g)
  end
  
  return res
end




