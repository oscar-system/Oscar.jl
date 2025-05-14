################################################################################
# Methods for generating Drinfeld-Hecke forms with a local strategy over a field of characteristic 0
#
# Over a field of characteristic 0, we have:
#
# For 1 != g ∈ G, there exists a Drinfeld-Hecke form with κ_g != 0 if and only if
#   (a) ker κ_g = V^g where V^g = ker(1 - M)
#   (b) codim(V^g) = 2
#   (c) det h⊥ = 1 for all h ∈ ZG(g) where h⊥ is h restricted to (V^g)⊥ = im(1 - M)
#
# For 1 ∈ G, the Drinfeld-Hecke forms are just the G-invariant bilinear alternating forms, i.e κ_1 with
#   (d) κ_1(v,w) = κ_1(gv,gw) for all g ∈ G and v,w ∈ V
# (Theorem 1.9 in Ram & Shepler: "Classification of graded Hecke algebras for complex reflection groups", 2002)
################################################################################

################################################################################
# Returns a dictionary mapping elements of G to matrices representing alternating bilinear forms
################################################################################
function generate_generic_forms_locally(G::MatrixGroup{T}, R::Field) where {T <: FieldElem}
  nonzero_forms = Dict{Tuple{MatrixGroupElem{T}, GroupConjClass}, MatElem}()
  number_of_parameters = 0
  
  # Iterate over the conjugacy classes and generate generic forms for the representatives
  for C in conjugacy_classes(G)
    g = representative(C)
    
    if is_one(g)
      κ_g = calculate_group_invariant_form(G, R)
    else
      κ_g = calculate_form_for_non_trivial_element(g, R)
    end
    
    # if κ_g != 0, we save it and add the number of parameters
    if !is_zero(κ_g)
      nonzero_forms[(g, C)] = κ_g
      number_of_parameters = number_of_parameters + ngens(base_ring(κ_g))
    end
  end

  # Initialize result
  forms = Dict{MatrixGroupElem{T}, MatElem}()

  # If there weren't any nonzero forms, return empty dictionary
  if length(nonzero_forms) == 0 return forms end

  # Now we know the number of parameters and can create our new base ring
  parameters = number_of_parameters == 1 ? ["t"] : ["t" * string(i) for i in 1:number_of_parameters]
  S, _ = polynomial_ring(R, parameters)
  
  # Next we shift all forms into S and calculate all remaining forms for the according class
  current_parameter_index = 1
  for ((g, C), κ_g) in classes_with_nonzero_forms
    # Shift κ_g into S using an homomorphism
    S_g = base_ring(κ_g)
    next_index = current_parameter_index + ngens(S_g)
    φ_g = hom(S_g, S, [S[i] for i in current_parameter_index:(next_index - 1)])
    κ_g = map(x -> φ_g(x), κ_g)
    current_parameter_index = next_index
    
    # Calculate all forms in conjugacy class and set everything to forms dictionary
    for c in C
      forms[c] = calculate_form_for_conjugate(g, c, κ_g)
    end
  end

  return forms
end

################################################################################
# Returns alternating bilinear form satisfying (a), (b), (c) represented as a matrix in the standard basis
################################################################################
function calculate_form_for_non_trivial_element(g::MatrixGroupElem{T}, R::Field) where {T <: FieldElem}
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

################################################################################
# Returns G-invariant alternating bilinear form represented as a matrix in the standard basis
################################################################################
function calculate_group_invariant_form(G::MatrixGroup{T}, R::Field) where {T <: FieldElem}
  M = build_group_invariant_relation_matrix(G)
  sol = solve_and_parametrize(M, R)
  
  # Convert solution to matrix
  result = zero_matrix(parent(sol[1]), n, n)
  
  for ((i,j),k) in map
    result[i,j] = sol[k]
    result[j,i] = -sol[k]
  end

  return result
end

################################################################################
# Translates the relations
#   κ(gu, gv) = κ(u, v) for all g ∈ G and u, v ∈ V
# into a matrix M representing the LES Mx = 0 such that x represents an alternating bilinear form
#
# For this note that
# - we only need to solve the relations for basis elements {v1,...,vn} of V
# - due to skew-symmetry it is enough to solve the relations for all combinations of i < j < k
# - if A = (a_ij) is the matrix corresponding to g in the given basis, then gvi = sum_l (a_li * vl). 
#
# With this the relations translate to
#   sum_{l < k} (a_li a_kj − a_ki a_lj) κ(vl,vk) − κ(vi,vj) = 0
# for all i < j
################################################################################
function build_group_invariant_relation_matrix(G::MatrixGroup)
  K = base_ring(G)
  n = degree(G)
  map = build_local_map(n)
  m = length(map)
  
  # Start to collect rows for relations
  rows = []
  
  for g in G
    # For g = 1 the relations are always true
    if is_one(g) continue end
  
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
  
  return (M, map)
end

################################################################################
# Calculate form for conjugate c = h−1gh of g by formula 
#   κ_h−1gh(v, w) = κ_g(hv, hw)
################################################################################
function calculate_form_for_conjugate(
  g::MatrixGroupElem{T},
  c::MatrixGroupElem{T},
  κ_g::MatElem
) where {T <: FieldElem}
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
# Generate index map that maps a pair (i,j) with i < j and 0 < i, j <= n to a fixed index
# Can be used for single group element to map solution vector to matrix and vice versa
################################################################################
function build_local_map(n::Int)
  map = Dict{Tuple{Int, Int}, Int}()
  
  k = 1
  for i in 1:n, j in (i+1):n
    map[(i,j)] = k
    k = k + 1
  end

  return map
end

