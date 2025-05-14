################################################################################
# Methods for generating Drinfeld-Hecke forms globally (all at once) 
#
# Let (κ_g)_g∈G be a family of alternating bilinear forms. Then
#   κ = sum_g∈G (κ_g * g)
# defines a Drinfeld-Hecke form if and only if
#   (a) κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all g ∈ G and u, v, w ∈ V
#   (a) κ_g(hu, hv) = κ_h^-1gh(u, v) for all g, h ∈ G and u, v ∈ V
# (Lemma 1.5 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
#
# These relations can be translated into a matrix M defining an LES of the form Mx = 0.
################################################################################

################################################################################
# Returns a parametrized family (κ_g)_g∈G of alternating bilinear forms
################################################################################
function generate_generic_forms_globally(G::MatrixGroup{T}, R::Ring) where {T <: FieldElem}
  M, map = build_relation_matrix(G)
  
  # Calculate global solution represented by a vector
  sol = solve_and_parametrize(M, R)
  
  # Convert solution to dictionary of matrices representing alternating bilinear forms
  S = parent(sol[1])
  m = length(sol)
  n = degree(G)
  
  # Initialize forms
  forms = Dict{MatrixGroupElem{T}, MatElem}()
  for g in G
    forms[g] = zero_matrix(S, n, n)
  end
  
  # Populate forms
  for ((g,i,j),k) in map
    forms[g][i,j] = sol[k]
    forms[g][j,i] = -sol[k]
  end

  return forms
end

