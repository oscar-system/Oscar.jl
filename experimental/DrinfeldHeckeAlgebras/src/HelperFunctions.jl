################################################################################
# Helper functions for generic Drinfeld-Hecke form generation and validation
################################################################################

################################################################################
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
################################################################################
function build_relation_matrix(G::MatrixGroup)
  K = base_ring(G)
  n = degree(G)
  map = build_global_map(G)
  m = length(map)
  
  # Start to collect rows for relations
  rows = []
  
 for g in G
    # First we add the relations (I) - (IV) for g != 1, since for g = 1, they are always true
    if !is_one(g)
      A = matrix(g)

      for i in 1:n, j in (i+1):n, k in (j+1):n
        for l in 1:n
          row = fill(K(), m)

          if l == k   # (II)
            row[map[(g,i,j)]] = A[k,k] - K(1)
            row[map[(g,j,k)]] = A[k,i]
            row[map[(g,i,k)]] = -A[k,j]

          elseif l == i   # (III)
            row[map[(g,i,j)]] = A[i,k]
            row[map[(g,j,k)]] = A[i,i] - K(1)
            row[map[(g,i,k)]] = -A[i,j]

          elseif l == j   # (IV)
            row[map[(g,i,j)]] = A[j,k]
            row[map[(g,j,k)]] = A[j,i]
            row[map[(g,i,k)]] = -A[j,j] + K(1)

          else  # (I)
            row[map[(g,i,j)]] = A[l,k]
            row[map[(g,j,k)]] = A[l,i]
            row[map[(g,i,k)]] = -A[l,j]
          end

          # If the row is nonzero, we add it to the LES
          if !is_zero(row) push!(rows, row) end
        end
      end
    end

    # Now we add the relations (V) for all h ∈ G and i < j
    for h in G
      # If h is one, they are always true
      if is_one(h) continue end

      c = inv(h) * g * h
      A = matrix(h)
      row = fill(K(), m)

      for i in 1:n, j in (i+1):n
        for l in 1:n, k in (l+1):n
          if g == c
            if l == i && k == j
              row[map[(g,l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j] - K(1)
            else
              row[map[(g,l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j]
            end
          else
            row[map[(g,l,k)]] = A[l,i] * A[k,j] - A[k,i] * A[l,j]

            if l == i && k == j
              row[map[(c,i,j)]] = K(-1)
            end
          end
        end
      end

      if !is_zero(row)
        push!(rows, row)
      end
    end
  end

  # Combine collected rows to matrix and return
  M = matrix(K, length(rows), m, vcat(rows...))
  
  return (M, map)
end

################################################################################
# Map elements g ∈ G and matrix indices i < j to global solution column index k
################################################################################
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

################################################################################
# Solves the LES Mx = 0 and returns a parametrized solution over the ring R
################################################################################
function solve_and_parametrize(M::MatElem{T}, R::Ring) where {T <: FieldElem}
  nullity, kernel_basis = nullspace(M)
  m = nrows(kernel_basis)
  
  # If nullity = 0, there is only the solution 0
  if nullity == 0
    return fill(R(), m)
  end
  
  # For creating a parametrized solution, we work over the polynomial ring S = R[t1,...tn]
  parameters = nullity == 1 ? ["t"] : ["t" * string(i) for i in 1:nullity]
  S, _ = polynomial_ring(R, parameters)

  # Use kernel basis to create a parametrized solution for a Drinfeld-Hecke form
  sol = fill(S(), m)
  
  for j in 1:nullity
    sol = sol + [S[j] * kernel_basis[i, j] for i in 1:m]
  end

  return sol
end
