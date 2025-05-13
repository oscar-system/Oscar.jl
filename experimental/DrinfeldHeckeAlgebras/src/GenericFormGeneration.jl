################################################################################
# Methods for generating and validating Drinfeld-Hecke forms
#
# General theory:
#
# Let (κ_g)_g∈G be a family of alternating bilinear forms. Then
#   κ = sum_g∈G (κ_g * g)
# defines a Drinfeld-Hecke form if and only if
#   (I)  κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all g ∈ G and u, v, w ∈ V
#   (II) κ_g(hu, hv) = κ_h^-1gh(u, v) for all g, h ∈ G and u, v ∈ V
# (Lemma 1.5 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
#
# Finding a solution for (I) for a representative of a conjugacy class, we can use (II) to find the remaining forms.
#
# Over a field of characteristic 0, instead of solving (I) we can also use:
#
# For 1 != g ∈ G, there exists a Drinfeld-Hecke form with κ_g != 0 if and only if
# (a) ker κ_g = V^g where V^g = ker(1 - M)
# (b) codim(V^g) = 2
# (c) det h⊥ = 1 for all h ∈ ZG(g) where h⊥ is h restricted to (V^g)⊥ = im(1 - M)
#
# For 1 ∈ G, the Drinfeld-Hecke forms are just the G-invariant bilinear alternating forms, i.e κ_1 with
#   κ_1(v,w) = κ_1(gv,gw) for all g ∈ G and v,w ∈ V
# (Theorem 1.9 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
################################################################################

################################################################################
# Drinfeld-Hecke form validation
################################################################################

# Check whether a set of bilinear forms defines a valid Drinfeld-Hecke form
function are_valid_forms(forms::Dict{MatrixGroupElem{T}, MatElem{S}}) where {T <: FieldElem, S <: RingElem}
  # If the forms are empty, they define the trivial (zero) Drinfeld-Hecke form
  if length(forms) == 0 return true end
  
  G = parent(first(forms)[1])
  n = degree(G)
  
  for C in conjugacy_classes(G)
    g = representative(C)
    
    if !haskey(forms, g)
      # If the form for g is zero, all other forms in this class must also be zero
      for h in C
        if haskey(forms, h) return false end
      end
    
    else
      # Otherwise we first need to check that the form for g defines a valid form
      κ_g = forms[g]
      if !is_valid_form(g, κ_g) return false end
      
      # Now we know that the form for g is valid. We need to check the other forms in this class
      for c in C
        # If κ_g is nonzero, then also κ_c needs to be nonzero
        if !haskey(forms, c) return false end
        
        # Check if conjugate form is valid
        if !is_valid_conjugate_form(g, κ_g, c, forms[c]) return false end
      end
    end
  end
  
  return true
end

# Check whether a single bilinear form fulfills the relations
#   κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all u, v, w ∈ V
function is_valid_form(g::MatrixGroupElem{T}, κ_g::MatElem{S}) where {T <: FieldElem, S <: RingElem}
  M = build_relation_matrix(g)
  sol = matrix_to_solution(κ_g)
  display(M)
  display(κ_g)
  display(M * sol)
  return is_zero(M * sol)
end

# Check whether κ_c fulfills the equation
#   κ_c(v, w) = κ_g(hv, hw)
# for a conjugate c = h−1gh of g
function is_valid_conjugate_form(
  g::MatrixGroupElem{T}, 
  c::MatrixGroupElem{T}, 
  κ_g::MatElem{S}, 
  κ_c::MatElem{S}
) where {T <: FieldElem, S <: RingElem}
  κ_c_correct = calculate_form_for_conjugate(g, c, κ_g)
  
  return κ_c_correct == κ_c
end

################################################################################
# Generic Drinfeld-Hecke form generation
################################################################################

function generate_generic_forms(G::MatrixGroup{T}, R::Ring) where {T <: FieldElem}
  # Iterate over the conjugacy classes and generate generic forms for the representatives
  classes = conjugacy_classes(G)
  classes_with_nonzero_forms = Dict{Tuple{MatrixGroupElem{T}, GroupConjClass}, MatElem}()
  number_of_parameters = 0

  for C in classes
    g = representative(C)
    
    # Generate generic form for representative
    κ_g = generate_generic_form(g, R)
    
    # if κ_g != 0, we save it and add the number of parameters
    if !is_zero(κ_g)
      classes_with_nonzero_forms[(g, C)] = κ_g
      number_of_parameters = number_of_parameters + ngens(base_ring(κ_g))
    end
  end

  # If there weren't any nonzero forms, return empty dictionary
  if length(classes_with_nonzero_forms) == 0
    return Dict{MatrixGroupElem{T}, MatElem{S}}()
  end

  # Now we know the number of parameters and can create our new base ring
  parameters = number_of_parameters == 1 ? ["t"] : ["t" * string(i) for i in 1:number_of_parameters]
  S, _ = polynomial_ring(R, parameters)
  forms = Dict{MatrixGroupElem{T}, MatElem{elem_type(typeof(S))}}()
  
  # Next we shift all forms into S and calculate all remaining forms for the according class
  current_parameter_index = 1
  for ((g, C), κ_g) in classes_with_nonzero_forms
    # Shift κ_g into S using an homomorphism
    S_g = base_ring(κ_g)
    next_index = current_parameter_index + ngens(S_g)
    φ_g = hom(S_g, S, [S[i] for i in current_parameter_index:(next_index - 1)])
    κ_g = map(x -> φ_g(x), κ_g)
    current_parameter_index = next_index
    
    # Calculate remaining forms in conjugacy class and set everything to forms dictionary
    for c in C
      forms[c] = calculate_form_for_conjugate(g, c, κ_g)
    end
  end

  return forms
end

################################################################################
# Conjugation handling
################################################################################

# Calculate form for conjugate c = h−1gh of g by formula 
#   κ_h−1gh(v, w) = κ_g(hv, hw)
function calculate_form_for_conjugate(
  g::MatrixGroupElem{T},
  c::MatrixGroupElem{T},
  κ_g::MatElem{S}
) where {T <: FieldElem, S <: RingElem}
  if c == g return κ_g end

  is_conj, h = is_conjugate_with_data(parent(g), g, c)
  
  if !is_conj
    throw(ArgumentError("Input c needs to be a conjugate of input g."))
  end
  
  A = matrix(h)
  n = nrows(g)
  κ_c = zero_matrix(base_ring(κ_g), n, n)
  
  for i in 1:n, j in (i+1):n
    # the i-th column of A is the result of letting h act on the i-th standard basis vector
    A_i = matrix(A[:,i])
    A_j = matrix(A[:,j])

    κ_c[i,j] = (transpose(A_i) * κ_g * A_j)[1,1]
    κ_c[j,i] = -κ_c[i,j]
  end

  return κ_c
end

################################################################################
# Generic Drinfeld-Hecke form generation for a group element
################################################################################

function generate_generic_form(g::MatrixGroupElem{T}, R::Ring) where {T <: FieldElem}
  if R isa Field && characteristic(R) == 0
    κ_g = generate_by_theorem_1_9(g, R)
  else
    κ_g = generate_by_lemma_1_5(g, R)
  end
end

################################################################################
# Generic Drinfeld-Hecke form generation according to Theorem 1.9
################################################################################

function generate_by_theorem_1_9(g::MatrixGroupElem{T}, R::Ring) where {T <: FieldElem}
  if is_one(g)
    return calculate_group_invariant_form(parent(g), R)
  else
    return calculate_form_for_non_trivial_element(g, R)
  end
end

# For 1 != g ∈ G, there exists a Drinfeld-Hecke form with κ_g != 0 if and only if
# (a) ker κ_g = V^g where V^g = ker(1 - M)
# (b) codim(V^g) = 2
# (c) det h⊥ = 1 for all h ∈ ZG(g) where h⊥ is h restricted to (V^g)⊥ = im(1 - M)
function calculate_form_for_non_trivial_element(g::MatrixGroupElem{T}, R::Ring) where {T <: FieldElem}
  G = parent(g)
  n = degree(G)
  
  I = identity_matrix(R, n)
  
  # Calculate dim and basis of V^g = ker(id - g)
  dim_Vg, basis_Vg = nullspace(I - matrix(g))
  
  # Check codim(V^g) = 2
  if n - dim_Vg != 2 
    return zero_matrix(R, n, n)  
  end
  
  # Find a basis for (V^g)⊥ = im(id - g)
  # Step 1: Extend basis_Vg to basis of V
  full_basis = basis_Vg
  for i in 1:n
    extended_basis = hcat(full_basis, matrix(I[:, i]))
    if rank(extended_basis) > rank(full_basis)
      full_basis = extended_basis
    end
  end
  
  # Step 2: Remove V^g basis to get complement basis
  complement_basis = full_basis[:, (dim_Vg + 1):n]
  
  # Apply I - matrix(g) to get a basis of (V^g)⊥ = im(id - g)
  basis_Vg⊥ = (I - matrix(g)) * complement_basis
  
  # Check det(h⊥) = 1 for elements in the centralizer restricted to (V^g)⊥
  ZG, _ = centralizer(G,g)
  
  # Compute left-inverse of basis_Vg⊥ (via Moore-Penrose pseudoinverse formula)
  left_inverse = inv(transpose(basis_Vg⊥) * basis_Vg⊥) * transpose(basis_Vg⊥)
  
  for h in ZG
    # Restrict h to the space (V^g)⊥
    h⊥ = left_inverse * matrix(h) * basis_Vg⊥
    
    if det(h⊥) != one(R) 
      return zero_matrix(R, n, n)  
    end
  end

  # We know that κ_g is now defined by its value on the basis {v1,v2} of (V^g)⊥, so we need one parameter
  S, _ = polynomial_ring(R, ["t"])
  result = zero_matrix(S, n, n)
  
  # Create matrix of κ_g in the basis {v1,v2}
  result[1,2] = S[1]
  result[2,1] = -S[1]
  
  # Base change to the standard basis
  combined_basis = hcat(basis_Vg, basis_Vg⊥)
  result = transpose(combined_basis) * result * combined_basis
  
  return result
end

# Calculate generic G-invariant bilinear alternating form, i.e κ with
#   κ(v,w) = κ(gv,gw) for all g ∈ G and v,w ∈ V
function calculate_group_invariant_form(G::MatrixGroup{T}, R::Ring) where {T <: FieldElem}
  K = base_ring(G)
  n = degree(G)
  
  map = build_map(n)
  m = length(map)
  rows = []
  
  for g in G
    # κ(gvi, gvj) = κ(vi, vj) is true if g = 1, so we can skip that case
    if is_one(g)
      continue
    end
  
    A = matrix(g)
    row = fill(K(), m)

    # The relations translate to 
    # sum_{l < k} (a_li a_kj − a_ki a_lj) κ(vl,vk) − κ(vi,vj) = 0
    for i in 1:n, j in (i+1):n
      for l in 1:n, k in (l+1):n
        if l == i && k == j
          row[map[(l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j] - K(1)
        else
          row[map[(l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j]
        end
      end
    end
  
    if !is_zero(row)
      push!(rows, row)
    end
  end

  # Combine collected rows to matrix
  M = matrix(K, length(rows), m, vcat(rows...))

  return solve_and_parametrize(M, R)
end

################################################################################
# Relation solving according to Lemma 1.5
################################################################################

# Translates for g ∈ G the relations
#   κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all u, v, w ∈ V
# into a matrix M representing the LES Mx = 0 and returns a solution
#
# For this note that
# - we only need to solve the relations for basis elements {v1,...,vn} of V
# - due to skew-symmetry it is enough to solve 
#   κ_g(vi, vj)(gvk − vk) + κ_g(vj, vk)(gvi − vi) - κ_g(vi, vk)(gvj − vj) = 0 
#   for all combinations of i < j < k
#
# If A = (a_ij) is the matrix corresponding to g in the given basis, then gvi = sum_l (a_li * vl). 
# Using this and rewriting the relations in the basis gives us the following relations for the coefficients:
#   (I)   κg(vi, vj) a_lk       + κg(vj, vk) a_li       - κg(vi, vk) a_lj      = 0    if l != i,j,k
#   (II)  κg(vi, vj) (a_kk − 1) + κg(vj, vk) a_ki       - κg(vi, vk) a_kj      = 0    (l == k)
#   (III) κg(vi, vj) a_ik       + κg(vj, vk) (a_ii − 1) - κg(vi, vk) a_ij      = 0    (l == i)
#   (IV)  κg(vi, vj) a_jk       + κg(vj, vk) a_ji       - κg(vi, vk) (ajj − 1) = 0    (l == j)
function generate_by_lemma_1_5(g::MatrixGroupElem{T}, R::Ring) where {T <: FieldElem}
  M = build_relation_matrix(g)

  return solve_and_parametrize(M, R)
end

# Translates for g ∈ G the relations
#   κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all u, v, w ∈ V
# into a matrix M representing the LES Mx = 0 such that x represents a solution for κ_g
#
# For this note that
# - we only need to solve the relations for basis elements {v1,...,vn} of V
# - due to skew-symmetry it is enough to solve 
#   κ_g(vi, vj)(gvk − vk) + κ_g(vj, vk)(gvi − vi) - κ_g(vi, vk)(gvj − vj) = 0 
#   for all combinations of i < j < k
#
# If A = (a_ij) is the matrix corresponding to g in the given basis, then gvi = sum_l (a_li * vl). 
# Using this and rewriting the relations in the basis gives us the following relations for the coefficients:
#   (I)   κg(vi, vj) a_lk       + κg(vj, vk) a_li       - κg(vi, vk) a_lj      = 0    if l != i,j,k
#   (II)  κg(vi, vj) (a_kk − 1) + κg(vj, vk) a_ki       - κg(vi, vk) a_kj      = 0    (l == k)
#   (III) κg(vi, vj) a_ik       + κg(vj, vk) (a_ii − 1) - κg(vi, vk) a_ij      = 0    (l == i)
#   (IV)  κg(vi, vj) a_jk       + κg(vj, vk) a_ji       - κg(vi, vk) (ajj − 1) = 0    (l == j)
function build_relation_matrix(g::MatrixGroupElem{T}) where {T <: FieldElem}
  K = base_ring(g)
  n = nrows(g)
  map = build_map(n)
  m = length(map)
  
  # If g = 1, the relations are always true
  if is_one(g)
    return zero_matrix(K, 0, m)
  end

  rows = []
  A = matrix(g)

  for i in 1:n, j in (i+1):n, k in (j+1):n
    for l in 1:n
      row = fill(K(), m)
      
      if l == k   # (II)
        row[map[(i,j)]] = A[k,k] - K(1)
        row[map[(j,k)]] = A[k,i]
        row[map[(i,k)]] = -A[k,j]
        
      elseif l == i   # (III)
        row[map[(i,j)]] = A[i,k]
        row[map[(j,k)]] = A[i,i] - K(1)
        row[map[(i,k)]] = -A[i,j]
        
      elseif l == j   # (IV)
        row[map[(i,j)]] = A[j,k]
        row[map[(j,k)]] = A[j,i]
        row[map[(i,k)]] = -A[j,j] + K(1)
        
      else  # (I)
        row[map[(i,j)]] = A[l,k]
        row[map[(j,k)]] = A[l,i]
        row[map[(i,k)]] = -A[l,j]
      end
    
      if !is_zero(row)
        push!(rows, row)
      end
    end
  end

  # Combine collected rows to matrix
  return matrix(K, length(rows), m, vcat(rows...))
end

################################################################################
# Helper functions
################################################################################

# Solves the LES Mx = 0 and returns a parametrized solution
function solve_and_parametrize(M::MatElem{T}, R::Ring) where {T <: FieldElem}
  nullity, kernel_basis = nullspace(M)
  m = nrows(kernel_basis)
  
  # If nullity = 0, there is only the solution 0
  if nullity == 0
    return solution_to_matrix(fill(R(), m))
  end
  
  # For creating a parametrized solution, we work over the polynomial ring S = R[t1,...tn]
  S, _ = polynomial_ring(R, ["t" * string(i) for i in 1:nullity])

  # Use kernel basis to create a parametrized solution for a Drinfeld-Hecke form
  sol = fill(S(), m)
  
  for j in 1:nullity
    sol = sol + [S[j] * kernel_basis[i, j] for i in 1:m]
  end

  return solution_to_matrix(sol)
end

# Convert solution to a matrix representing an alternating bilinear form
function solution_to_matrix(sol::Vector)
  m = length(sol)
  n = Int((1 + sqrt(1 + 8*m)) / 2)
  
  result = zero_matrix(parent(sol[1]), n, n)
  
  for ((i,j),k) in build_map(n)
    result[i,j] = sol[k]
    result[j,i] = -sol[k]
  end

  return result
end

# Convert matrix representing an alternating bilinear form to a (potential solution) vector
function matrix_to_solution(mat::MatElem)
  map = build_map(nrows(mat))
  m = length(map)
  sol = fill(base_ring(mat)(), m)

  for ((i,j),k) in map
    sol[k] = mat[i,j]
  end

  return sol
end

# Generate index map that maps a pair (i,j) with i < j and 0 < i, j <= n to a fixed index
function build_map(n::Int)
  map = Dict{Tuple{Int, Int}, Int}()
  
  k = 1
  for i in 1:n, j in (i+1):n
    map[(i,j)] = k
    k = k + 1
  end

  return map
end

