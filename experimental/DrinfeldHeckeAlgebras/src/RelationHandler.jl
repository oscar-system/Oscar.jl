################################################################################
# Struct for relation handling of Drinfeld-Hecke forms
################################################################################

# TODO: Haben diese Elemente einen Namen? 1.9 Ram Shepler
# TODO: Gibt es in Characteristic 0 ein Beispiel bei dem es nur die 0-DH form gibt?
# TODO: Gibt es Drinfeld-Hecke algebren die keine Rational-Cherednik algebren sind für complexe Spiegelungsgruppe?
# TODO: Gibt es DH algebren die keine Symplectic reflection algebras sind?
# TODO. generic instead of parametrized

# Handles the relations below by translating them into a matrix M defining an LES of the form Mx = 0
# κg(u, v)(gw − w) + κg(v, w)(gu − u) + κg(w, u)(gv − v) = 0 for all g ∈ G and u, v, w ∈ V
# κg(hu, hv) = κh^-1gh(u, v) for all g, h ∈ G and u, v ∈ V
mutable struct RelationHandler{T <: RingElem}
  group::MatrixGroup{T}
  map::Dict{Tuple{MatrixGroupElem{T}, Int, Int}, Int}
  relation_matrix::MatrixElem{T}
  
  # For each g ∈ G, κg is defined by its values on the basis elements of V. 
  # So if v1, ... vn is a basis of V, we need to calculate κg(vi,vj) for all combinations of i and j
  # Since all κg are alternating bilinear forms, it is enough to look at i < j
  # We will set up an LES of the form Mx = 0 where each column of M corresponds to a value κg(vi,vj)
  function RelationHandler(G::MatrixGroup{T}) where {T <: RingElem}
    handler = new{T}()
    handler.group = G
    
    d = degree(G)
    
    # Generate index map that maps a triple (g,i,j) with i < j to a column index
    handler.map = Dict{Tuple{MatrixGroupElem{T}, Int, Int}, Int}()
    
    k = 1
    for g in G, i in 1:d, j in (i+1):d
      handler.map[(g,i,j)] = k
      k = k + 1
    end
    
    # Translates the relations to a matrix M defining an LES Mx = 0
    R = base_ring(G)
    m = length(handler.map)
    
    # Start collecting the rows of M
    rows = []
    conjugates = Dict()
    
    for g in G
      A = matrix(g)
      
      # We add the relations κg(vi, vj)(gvk − vk) + κg(vj, vk)(gvi − vi) - κg(vi, vk)(gvj − vj) = 0
      # for all combinations of i < j < k
      for i in 1:d, j in (i+1):d, k in (j+1):d
        # The relations are true if g = 1, so we can skip that case
        if is_one(g)
          continue
        end
        
        # Note that if A = (a_ij) is the matrix corresponding to g in the given basis, then gvi = sum_l (a_li * vl). 
        # Using this and rewriting the relations in the basis gives us the following relations for the coefficients:
        #   (I)   κg(vi, vj) a_lk       + κg(vj, vk) a_li       - κg(vi, vk) a_lj      = 0    if l != i,j,k
        #   (II)  κg(vi, vj) (a_kk − 1) + κg(vj, vk) a_ki       - κg(vi, vk) a_kj      = 0    (l == k)
        #   (III) κg(vi, vj) a_ik       + κg(vj, vk) (a_ii − 1) - κg(vi, vk) a_ij      = 0    (l == i)
        #   (IV)  κg(vi, vj) a_jk       + κg(vj, vk) a_ji       - κg(vi, vk) (ajj − 1) = 0    (l == j)
        for l in 1:d
          if l == k   # (II)
            row = fill(R(), m)
            
            row[handler.map[(g,i,j)]] = A[k,k] - R(1)
            row[handler.map[(g,j,k)]] = A[k,i]
            row[handler.map[(g,i,k)]] = -A[k,j]
            
          elseif l == i   # (III)
            row = fill(R(), m)
            
            row[handler.map[(g,i,j)]] = A[i,k]
            row[handler.map[(g,j,k)]] = A[i,i] - R(1)
            row[handler.map[(g,i,k)]] = -A[i,j]
            
          elseif l == j   # (IV)
            row = fill(R(), m)
            
            row[handler.map[(g,i,j)]] = A[j,k]
            row[handler.map[(g,j,k)]] = A[j,i]
            row[handler.map[(g,i,k)]] = -A[j,j] + R(1)
            
          else  # (I)
            row = fill(R(), m)
            
            row[handler.map[(g,i,j)]] = A[l,k]
            row[handler.map[(g,j,k)]] = A[l,i]
            row[handler.map[(g,i,k)]] = -A[l,j]
          end
        
          if !is_zero(row)
            push!(rows, row)
          end
        end
      end
    
      # We add the relations κg(hvi, hvj) = κh^-1gh(vi, vj) for all h ∈ G and i < j
      for h in G, i in 1:d, j in (i+1):d
        # The relations are true if h = 1, so we can skip that case
        if is_one(h)
          continue
        end
      
        B = matrix(h)

        if haskey(conjugates, (g,h))
          c = conjugates[(g,h)]
        else
          c = conj(g, h)
          conjugates[(g,h)] = c
        end
      
        row = fill(R(), m)
        
        # The relations translate to 
        # sum_{l < k} (bli bkj − bki blj) κg(vl,vk) − κh−1gh(vi,vj) = 0
        # where B = (b_ij) is the matrix corresponding to h
        for l in 1:d, k in (l+1):d
          if g == c && l == i && k == j
            row[handler.map[(g,l,k)]] = B[l,i] * B[k,j] - B[k,i] * B[l,j] - R(1)
          else
            row[handler.map[(g,l,k)]] = B[l,i] * B[k,j] - B[k,i] * B[l,j]
            row[handler.map[(c,i,j)]] = R(-1)
          end
          # TODO I think this is not correct yet
        end
      
        if !is_zero(row)
          push!(rows, row)
        end
      end
    end
  
    # Combine collected rows to matrix
    handler.relation_matrix = matrix(R, length(rows), m, vcat(rows...))
    
    return handler
  end
end

################################################################################
# Generic functions
################################################################################

relation_matrix(handler::RelationHandler) = handler.relation_matrix
map(handler::RelationHandler) = handler.map
group(handler::RelationHandler) = handler.group

# Over characteristic 0, for 1 != g ∈ G, there exists a Drinfeld-Hecke form κg != 0 if and only if
# (a) ker κg = V^g where V^g = ker(1 - M)
# (b) codim(V^g) = 2
# (c) det h⊥ = 1 for all h ∈ ZG(g) where h⊥ is h restricted to (V^g)⊥ = im(1 - M)
#  
# For 1 ∈ G, the Drinfeld-Hecke forms are just the G-invariant bilinear alternating forms, i.e κ1 with
# κ1(v,w) = κ1(gv,gw) for all g ∈ G and v,w ∈ V
# (Theorem 1.9, 2002 Ram & Shepler: Classification of graded Hecke algebras for complex reflection groups)
function generate_forms_for_conjugacy_classes(G::MatrixGroup{T}) where {T <: RingElem}
  R = base_ring(G)
  n = degree(G)
  I = identity_matrix(R, n)
  
  classes = conjugacy_classes(G)
  
  # keep track of conjugacy classes with nonzero forms and attach data
  # (dim_V^g, basis_V^g, dim_(V^g)⊥, basis_(V^g)⊥)
  classes_with_nonzero_forms = Dict{GroupConjClass, Tuple{Int, MatrixElem{T}, Int, MatrixElem{T}}}()
    
  for C in classes
    g = representative(C)
    
    # continue if g = 1
    if is_one(g) continue end
    
    # Calculate dim and basis of V^g = ker(id - g)
    dim_Vg, basis_Vg = nullspace(I - matrix(g))
    
    # Check dim((V^g)⊥) = codim(V^g) = 2
    dim_Vg⊥ = n - dim_Vg
    if dim_Vg⊥ != 2 continue end
    
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
      
      if det(h⊥) != one(R) continue end
    end
    
    # Now we know that (b) and (c) are satisfied.
    # Hence for this conjugacy class there exist nonzero forms corresponding to a Drinfeld-Hecke form
    classes_with_nonzero_forms[C] = (dim_Vg, basis_Vg, dim_Vg⊥, basis_Vg⊥)
  end

  d = length(classes_with_nonzero_forms)
  
  # Next we need to solve κ1(gvi, gvj) = κ1(vi, vj) for all g ∈ G and i < j to get all κ1
  # Generate index map that maps a pair (i,j) with i < j to a column index
  map = Dict{Tuple{Int, Int}, Int}()
  
  k = 1
  for i in 1:n, j in (i+1):n
    map[(i,j)] = k
    k = k + 1
  end

  m = length(map)
  rows = []
  
  for g in G
    # κ1(gvi, gvj) = κ1(vi, vj) is true if g = 1, so we can skip that case
    if is_one(g)
      continue
    end
  
    B = matrix(g)
    row = fill(R(), m)
    
    # The relations translate to 
    # sum_{l < k} (bli bkj − bki blj) κ1(vl,vk) − κ1(vi,vj) = 0
    for i in 1:n, j in (i+1):n
      for l in 1:n, k in (l+1):n
        if l == i && k == j
          row[map[(l,k)]] = B[l,i] * B[k,j] - B[k,i] * B[l,j] - R(1)
        else
          row[map[(l,k)]] = B[l,i] * B[k,j] - B[k,i] * B[l,j]
        end
      end
    end
  
    if !is_zero(row)
      push!(rows, row)
    end
  end

  # Combine collected rows to matrix
  M = matrix(R, length(rows), m, vcat(rows...))

  # Solve the LES Mx = 0
  nullity, K = nullspace(M)

  # From the theorem we know that the number of parameters is d + nullity
  dim = d + nullity
  parameters = dim == 1 ? ["t"] : ["t" * string(i) for i in 1:dim]
  S, _ = polynomial_ring(R, parameters)

  # We start to build the parametrized Drinfeld-Hecke forms
  forms = Dict{MatrixGroupElem, MatElem}()
  
  # First for 1 ∈ G
  if nullity > 0
    x = fill(S(), m)
    for j in 1:nullity
      x = x + [S[j] * K[i, j] for i in 1:m]
    end

    # Convert to a matrix
    B = zero_matrix(S, n, n)
    for ((i,j),k) in map
      B[i,j] = x[k]
      B[j,i] = -x[k]
    end
  
    forms[one(G)] = B
  end

  # Now for the conjugacy classes
  parameter_index = nullity + 1
  for (C, (dim_Vg, basis_Vg, dim_Vg⊥, basis_Vg⊥)) in classes_with_nonzero_forms
    g = representative(C)
    
    # Create matrix of κ_g in the basis (basis_Vg, basis_Vg⊥)
    B = zero_matrix(S, n, n)
    B[1,2] = S[parameter_index]
    B[2,1] = -S[parameter_index]

    # Base change to the standard basis
    combined_basis = hcat(basis_Vg, basis_Vg⊥)
    forms[g] = transpose(combined_basis) * B * combined_basis
    
    # Now calculate all remaining forms for the conjugacy class
    for h in C
      if h == g continue end
      
      H = matrix(h)
      D = zero_matrix(S, n, n)
      for i in 1:n, j in (i+1):n
        h_ei = H * I[:,i]
        h_ej = H * I[:,j]
        D[i,j] = (transpose(matrix(h_ei)) * forms[g] * matrix(h_ej))[1,1]
        D[j,i] = -D[i,j]
      end
    
      forms[h] = D
    end
    
    parameter_index = parameter_index + 1
  end
    
  return forms
end

