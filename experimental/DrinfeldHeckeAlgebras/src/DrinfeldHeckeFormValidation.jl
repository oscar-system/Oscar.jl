################################################################################
# Methods for validating Drinfeld-Hecke forms globally (all at once) 
#
# We check that for a family (κ_g)_g∈G of alternating bilinear forms the relations
#   (a) κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all g ∈ G and u, v, w ∈ V
#   (a) κ_g(hu, hv) = κ_h^-1gh(u, v) for all g, h ∈ G and u, v ∈ V
# are fulfilled. If they are, κ = sum_g∈G (κ_g * g) is a Drinfeld-Hecke form
#
# (Lemma 1.5 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
################################################################################

################################################################################
# Returns true if the given forms define a Drinfeld-Hecke algebra, false otherwise
################################################################################
function is_drinfeld_hecke_form(forms::Dict{MatrixGroupElem{T}, MatElem{S}}) where {T <: FieldElem, S <: RingElem}
  # If the forms are empty, they define the trivial (zero) Drinfeld-Hecke form
  if length(forms) == 0 return true end
  
  # Otherwise we extract group and ring data
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
