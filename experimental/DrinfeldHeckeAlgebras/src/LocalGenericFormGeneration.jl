#######################################
# Methods for generating Drinfeld-Hecke forms with a local strategy over a field of characteristic 0
#
# Over a field of characteristic 0, we have:
#
# For 1 != g ∈ G, there exists a Drinfeld-Hecke form with kappa_g != 0 if and only if
#   (a) ker kappa_g = V^g where V^g = ker(1 - M)
#   (b) codim(V^g) = 2
#   (c) det h⊥ = 1 for all h ∈ ZG(g) where h⊥ is h restricted to (V^g)⊥ = im(1 - M)
#
# For 1 ∈ G, the Drinfeld-Hecke forms are just the G-invariant bilinear alternating forms, i.e kappa_1 with
#   (d) kappa_1(v,w) = kappa_1(gv,gw) for all g ∈ G and v,w ∈ V
# (Theorem 1.9 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
#
# Cassandra Koenen, 2025
#######################################

#######################################
# Returns a dictionary mapping elements of G to matrices representing alternating bilinear forms
#######################################
function generate_generic_forms_locally(G::MatrixGroup{T}, R::Field) where {T <: FieldElem}
  nonzero_forms = Dict{Tuple{MatrixGroupElem{T}, GroupConjClass}, MatElem}()
  number_of_parameters = 0
  
  # Iterate over the conjugacy classes and generate generic forms for the representatives
  for C in conjugacy_classes(G)
    g = representative(C)
    
    if is_one(g)
      kappa_g = calculate_generic_group_invariant_form(G, R)
    else
      kappa_g = calculate_generic_form_for_non_trivial_element(g, R)
    end
    
    # if kappa_g != 0, we save it and add the number of parameters
    if !is_zero(kappa_g)
      nonzero_forms[(g, C)] = kappa_g
      number_of_parameters = number_of_parameters + ngens(base_ring(kappa_g))
    end
  end

  # If there weren't any nonzero forms, return empty dictionary
  if length(nonzero_forms) == 0 
    return Dict{MatrixGroupElem{T}, MatElem{elem_type(typeof(R))}}() 
  end

  # Now we know the number of parameters and can create our new base ring
  parameters = number_of_parameters == 1 ? ["t"] : ["t" * string(i) for i in 1:number_of_parameters]
  S, _ = polynomial_ring(R, parameters)
  forms = Dict{MatrixGroupElem{T}, MatElem{elem_type(typeof(S))}}()
  
  # Next we shift all forms into S and calculate all remaining forms for the according class
  current_parameter_index = 1
  for ((g, C), kappa_g) in nonzero_forms
    # Shift kappa_g into S using an homomorphism
    S_g = base_ring(kappa_g)
    next_index = current_parameter_index + ngens(S_g)
    phi_g = hom(S_g, S, [S[i] for i in current_parameter_index:(next_index - 1)])
    kappa_g = map(x -> phi_g(x), kappa_g)
    current_parameter_index = next_index
    
    # Calculate all forms in conjugacy class and set everything to forms dictionary
    for c in C
      forms[c] = calculate_form_for_conjugate(g, c, kappa_g)
    end
  end

  return forms
end

#######################################
# Returns alternating bilinear form satisfying (a), (b), (c) represented as a matrix in the standard basis
#######################################
function calculate_generic_form_for_non_trivial_element(g::MatrixGroupElem{T}, K::Field) where {T <: FieldElem}
  form = calculate_form_for_non_trivial_element(g, K)

  # We know that kappa_g is now defined by its value on the basis {v1,v2} of (V^g)⊥, so we need one parameter
  R, _ = polynomial_ring(K, ["t"])
  
  # Multiply form with t to parametrize it
  result = form * R[1]
  
  # Normalize for better readability
  n = degree(parent(g))
  for j in 1:n, i in 1:n
    if !is_zero(result[i,j])
      result = result / leading_coefficient(result[i,j])
    end
  end
  
  return result
end

#######################################
# Returns alternating bilinear form satisfying (a), (b), (c) represented as a matrix in the standard basis
#######################################
function calculate_form_for_non_trivial_element(g::MatrixGroupElem{T}, K::Field) where {T <: FieldElem}
  G = parent(g)
  n = degree(G)
  
  I = identity_matrix(K, n)
  A = I - matrix(g)
  
  # Calculate dim and basis of V^g = ker(id - g)
  dim_Vg, basis_Vg = nullspace(A)
  
  # Check codim(V^g) = 2
  if n - dim_Vg != 2 
    return zero_matrix(K, n, n)  
  end
  
  # We know that dim(V^g)⊥ = 2, hence we get a basis of (V^g)⊥ = im(id - g) by picking 2 columns of
  # id - g that are linearly independent 
  basis_Vg⊥ = I - matrix(g)
  
  for i in 1:n, j in (i+1):n
    basis_Vg⊥ = matrix(hcat(A[:,i], A[:,j]))

    # If the rank is 2, we stop
    if rank(basis_Vg⊥) == 2 break end
  end

  # Check det(h⊥) = 1 for elements in the centralizer restricted to (V^g)⊥
  ZGg, _ = centralizer(G,g)

  # Compute Moore-Penrose pseudoinverse of basis_Vg⊥
  basis_Vg⊥_pseudoinverse = inv(transpose(basis_Vg⊥) * basis_Vg⊥) * transpose(basis_Vg⊥)
  
  for h in ZGg
    # Restrict h to the space (V^g)⊥
    h⊥ = basis_Vg⊥_pseudoinverse * matrix(h) * basis_Vg⊥
    
    if det(h⊥) != one(K) 
      return zero_matrix(K, n, n)  
    end
  end

  # Create matrix of kappa_g in the basis {v1,v2}
  result = zero_matrix(K, n, n)
  result[1,2] = 1
  result[2,1] = -1
  
  # Base change to the standard basis
  B = hcat(basis_Vg⊥, basis_Vg)
  B_inv = inv(B)
  result = transpose(B_inv) * result * B_inv
  
  return result
end

#######################################
# Returns G-invariant alternating bilinear form represented as a matrix in the standard basis
#######################################
function calculate_generic_group_invariant_form(G::MatrixGroup{T}, R::Field) where {T <: FieldElem}
  M, map = build_group_invariant_relation_matrix(G)
  sol = solve_and_parametrize(M, R)
  n = degree(G)
  
  # Convert solution to matrix
  result = zero_matrix(parent(sol[1]), n, n)
  
  for ((i,j),k) in map
    result[i,j] = sol[k]
    result[j,i] = -sol[k]
  end

  return result
end

#######################################
# Translates the relations
#   kappa(gu, gv) = kappa(u, v) for all g ∈ G and u, v ∈ V
# into a matrix M representing the LES Mx = 0 such that x represents an alternating bilinear form
#
# For this note that
# - we only need to solve the relations for basis elements {v1,...,vn} of V
# - due to skew-symmetry it is enough to solve the relations for all combinations of i < j
# - if A = (a_ij) is the matrix corresponding to g in the given basis, then gvi = sum_l (a_li * vl). 
#
# With this the relations translate to
#   sum_{l < k} (a_li a_kj − a_ki a_lj) kappa(vl,vk) − kappa(vi,vj) = 0
# for all i < j
#######################################
function build_group_invariant_relation_matrix(G::MatrixGroup)
  K = base_ring(G)
  n = degree(G)
  map = build_local_map(n)
  m = length(map)
  
  # Start to collect equations for relations
  M = zero_matrix(K, 0, m)

  for g in gens(G)
    # For g = 1 the relations are always true
    if is_one(g) continue end
  
    A = matrix(g)
    row = zero_matrix(K, 1, m)
  
    # The relations translate to the equation
    # sum_{l < k} (a_li a_kj − a_ki a_lj) kappa(vl,vk) − kappa(vi,vj) = 0
    # for each i < j
    for i in 1:n, j in (i+1):n
      # Build the equation
      for l in 1:n, k in (l+1):n
        if l == i && k == j 
          row[1,map[(i,j)]] = A[i,i] * A[j,j] - A[j,i] * A[i,j] - K(1)
        else 
          row[1,map[(l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j]
        end
      end
    
      # Add the equation as a row
      if !is_zero(row)
        M = vcat(M, row)
      end
    end
  end

  return (M, map)
end

#######################################
# Calculate form for conjugate c = h−1gh of g by formula 
#   kappa_h−1gh(v, w) = kappa_g(hv, hw)
#######################################
function calculate_form_for_conjugate(
  g::MatrixGroupElem{T},
  c::MatrixGroupElem{T},
  kappa_g::MatElem
) where {T <: FieldElem}
  if c == g return kappa_g end

  is_conj, h = is_conjugate_with_data(parent(g), g, c)
  
  if !is_conj
    throw(ArgumentError("Input c needs to be a conjugate of input g."))
  end
  
  A = matrix(h)
  n = nrows(g)
  kappa_c = zero_matrix(base_ring(kappa_g), n, n)
  
  for i in 1:n, j in (i+1):n
    # the i-th column of A is the result of letting h act on the i-th standard basis vector
    A_i = matrix(A[:,i])
    A_j = matrix(A[:,j])

    kappa_c[i,j] = (transpose(A_i) * kappa_g * A_j)[1,1]
    kappa_c[j,i] = -kappa_c[i,j]
  end

  return kappa_c
end

#######################################
# Generate index map that maps a pair (i,j) with i < j and 0 < i, j <= n to a fixed index
# Can be used for single group element to map solution vector to matrix and vice versa
#######################################
function build_local_map(n::Int)
  map = Dict{Tuple{Int, Int}, Int}()
  
  k = 1
  for i in 1:n, j in (i+1):n
    map[(i,j)] = k
    k = k + 1
  end

  return map
end

