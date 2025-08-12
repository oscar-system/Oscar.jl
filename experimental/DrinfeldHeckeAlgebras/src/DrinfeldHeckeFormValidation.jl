#######################################
# Methods for validating Drinfeld-Hecke forms
#
# Cassandra Koenen, 2025
#######################################

#######################################
# Returns true if the given forms define a Drinfeld-Hecke algebra, false otherwise
#######################################
function is_drinfeld_hecke_form(forms::Dict{MatrixGroupElem{T}, MatElem{S}}) where {T <: FieldElem, S <: RingElem}
    # If the forms are empty, they define the trivial (zero) Drinfeld-Hecke form
    if length(forms) == 0 return true end
    
    # Otherwise we extract ring data
    g, κ_g = first(forms)
    R = base_ring(κ_g)
    
    if R isa Field && characteristic(R) == 0
      return is_drinfeld_hecke_form_local_strategy(forms)
    else
      return is_drinfeld_hecke_form_global_strategy(forms)
    end
end

#######################################
# Checks if the given forms define a valid Drinfeld-Hecke form using the local strategy 
# based upon the theorem referenced below
#
# (Theorem 1.9 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
#######################################
function is_drinfeld_hecke_form_local_strategy(
  forms::Dict{MatrixGroupElem{T}, MatElem{S}}
) where {T <: FieldElem, S <: RingElem}
  # Extract group and ring data
  g, κ_g = first(forms)
  G = parent(g)
  R = base_ring(κ_g)

  # We iterate over the conjugacy classes
  for C in conjugacy_classes(G)
    g = representative(C)
    
    # If g = 1, we only check G-invariancy
    if is_one(g)
      if haskey(forms, g) && !is_zero(forms[g])
        M, map = build_group_invariant_relation_matrix(G)
        
        # Translate 1-form into vector using the map
        sol = fill(R(), length(map))
        for ((i,j),k) in map
          sol[k] = forms[g][i,j]
        end
      
        # Check if sol is indeed a solution
        if !is_zero(M * sol)
          return false
        end
      end
    
      # Move on to the next class
      continue
    end

    # For nontrivial classes, check if there are nonzero forms in conjugacy class
    has_nonzero_form = false
    for h in C
      has_nonzero_form = has_nonzero_form || haskey(forms, h) && !is_zero(forms[h])
    end
  
    # If there aren't any nonzero forms, move on to the next class
    if !has_nonzero_form
      continue
    end
  
    # Otherwise check that there aren't any zero forms for this class
    for h in C
      if !haskey(forms, h) || is_zero(forms[h])
        return false
      end
    end
  
    # Now we calculate a valid form for g according to Ram & Shepler
    valid_form = calculate_form_for_non_trivial_element(g, R)
    
    # Then forms[g] must be a scalar multiple of valid_form
    n = degree(G)
    scalar = 1
    
    for i in 1:n, j in 1:n
      if !is_zero(valid_form[i,j])
        scalar = forms[g][i,j] / valid_form[i,j]
      end
    end
  
    if scalar * valid_form != forms[g]
      return false
    end
  
    # It remains to check the other elements in the conjugacy class
    for c in C
      if c == g
        continue
      end

      if forms[c] != calculate_form_for_conjugate(g, c, forms[g])
        return false
      end
    end
  end

  return true
end

#######################################
# Checks if the given forms define a valid Drinfeld-Hecke form using the global strategy
# of translating necessary relations into an LES 
#
# (Lemma 1.5 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
#######################################
function is_drinfeld_hecke_form_global_strategy(
  forms::Dict{MatrixGroupElem{T}, MatElem{S}}
) where {T <: FieldElem, S <: RingElem}
  # Extract group and ring data
  g, κ_g = first(forms)
  G = parent(g)
  R = base_ring(κ_g)
  
  # Build relation matrix
  M, map = build_relation_matrix(G)
  
  # Translate forms into vector using the map
  sol = fill(R(), length(map))
  for ((g,i,j),k) in map
    if haskey(forms, g)
      sol[k] = forms[g][i,j]
    end
  end

  # Check if sol is indeed a solution
  return is_zero(M * sol)
end
