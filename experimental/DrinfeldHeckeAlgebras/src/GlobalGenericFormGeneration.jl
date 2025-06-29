#######################################
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
#
# Cassandra Koenen, 2025
#######################################

#######################################
# Returns a parametrized family (κ_g)_g∈G of alternating bilinear forms defining a Drinfeld-Hecke form
#######################################
function generate_generic_forms_globally(G::MatrixGroup{T}, R::Ring) where {T <: FieldElem}
  M, map = build_relation_matrix(G)
  
  # Calculate global solution represented by a vector
  sol = solve_and_parametrize(M, R)
  
  # Convert solution to dictionary of matrices representing alternating bilinear forms
  S = parent(sol[1])
  m = length(sol)
  n = degree(G)
  
  # Initialize forms
  forms = Dict{MatrixGroupElem{T}, MatElem{elem_type(typeof(S))}}()
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

#######################################
# Translates the relations
#   (a) κ_g(u, v)(gw − w) + κ_g(v, w)(gu − u) + κ_g(w, u)(gv − v) = 0 for all g ∈ G and u, v, w ∈ V
#   (b) κ_g(hu, hv) = κ_h^-1gh(u, v) for all g, h ∈ G and u, v ∈ V
# into a matrix M representing the LES Mx = 0 such that x represents a global solution of all κ_g
#
# For this note that
# - we only need to solve the relations for basis elements {v1,...,vn} of V
# - due to skew-symmetry it is enough to solve the relations for all combinations of i < j < k
# - if A = (a_ij) is the matrix corresponding to g in the given basis, then gvi = sum_l (a_li * vl). 
#
# Using this we can rewrite the relations in (a) to
#   (I)   κ_g(vi, vj) a_lk       + κ_g(vj, vk) a_li       - κ_g(vi, vk) a_lj      = 0    if l != i,j,k
#   (II)  κ_g(vi, vj) (a_kk − 1) + κ_g(vj, vk) a_ki       - κ_g(vi, vk) a_kj      = 0    (l == k)
#   (III) κ_g(vi, vj) a_ik       + κ_g(vj, vk) (a_ii − 1) - κ_g(vi, vk) a_ij      = 0    (l == i)
#   (IV)  κ_g(vi, vj) a_jk       + κ_g(vj, vk) a_ji       - κ_g(vi, vk) (ajj − 1) = 0    (l == j)
# for all i < j < k and l and where A = (a_ij) is the matrix corresponding to g
#
# The relations in (b) on the other hand translate to
#   (V) sum_{l < k} (a_li a_kj − a_ki a_lj) κ_g(vl,vk) − κ_h−1gh(vi,vj) = 0
# for all i < j and where A = (a_ij) is the matrix corresponding to h
#######################################
function build_relation_matrix(G::MatrixGroup)
  K = base_ring(G)
  n = degree(G)
  map = build_global_map(G)
  m = length(map)
  
  # Start to collect equations for relations
  M = zero_matrix(K, 0, m)
  
  for g in G
    # First we add the relations (I) - (IV) for g != 1, since for g = 1, they are always true
    if !is_one(g)
      A = matrix(g)

      for i in 1:n, j in (i+1):n, k in (j+1):n
        for l in 1:n
          row = zero_matrix(K, 1, m)

          if l == k   # (II)
            row[1,map[(g,i,j)]] = A[k,k] - K(1)
            row[1,map[(g,j,k)]] = A[k,i]
            row[1,map[(g,i,k)]] = -A[k,j]

          elseif l == i   # (III)
            row[1,map[(g,i,j)]] = A[i,k]
            row[1,map[(g,j,k)]] = A[i,i] - K(1)
            row[1,map[(g,i,k)]] = -A[i,j]

          elseif l == j   # (IV)
            row[1,map[(g,i,j)]] = A[j,k]
            row[1,map[(g,j,k)]] = A[j,i]
            row[1,map[(g,i,k)]] = -A[j,j] + K(1)

          else  # (I)
            row[1,map[(g,i,j)]] = A[l,k]
            row[1,map[(g,j,k)]] = A[l,i]
            row[1,map[(g,i,k)]] = -A[l,j]
          end

          # If the row is nonzero, we add it to the LES
          if !is_zero(row) 
            M = vcat(M, row)  
          end
        end
      end
    end

    # Now we add the relations (V) for all h ∈ G and i < j
    for h in G
      # If h is one, they are always true
      if is_one(h) continue end

      c = inv(h) * g * h
      A = matrix(h)
      
      # We need
      # sum_{l < k} (a_li a_kj − a_ki a_lj) κ_g(vl,vk) − κ_h−1gh(vi,vj) = 0
      # for each i < j
      for i in 1:n, j in (i+1):n
        row = zero_matrix(K, 1, m)
        
        # Build the equation
        if g == c
          for l in 1:n, k in (l+1):n
            if l == i && k == j 
              row[1,map[(g,i,j)]] = A[i,i] * A[j,j] - A[j,i] * A[i,j] - K(1)
            else 
              row[1,map[(g,l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j]
            end
          end
        
        else
          row[1,map[(c,i,j)]] = -K(1)

          for l in 1:n, k in (l+1):n
            row[1,map[(g,l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j]
          end
        end
      
        if !is_zero(row)
          M = vcat(M, row)
        end
      end
    end
  end

  return (M, map)
end

#######################################
# Map elements g ∈ G and matrix indices i < j to global solution column index k
#######################################
function build_global_map(G::MatrixGroup)
  map = Dict{Tuple{MatrixGroupElem, Int, Int}, Int}()
  n = degree(G)
  
  k = 1
  for g in G, i in 1:n, j in (i+1):n
    map[(g,i,j)] = k
    k = k + 1
  end

  return map
end

