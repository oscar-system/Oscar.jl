################################################################################
# Struct for relation handling of Drinfeld-Hecke forms
################################################################################

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
            
            row[handler.map[(i,j)]] = A[l,k]
            row[handler.map[(j,k)]] = A[l,i]
            row[handler.map[(k,i)]] = -A[l,j]
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
        c = (h^-1) * g * h
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
