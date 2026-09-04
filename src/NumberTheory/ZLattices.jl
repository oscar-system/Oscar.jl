function _overlattice_orbits(L::ZZLat, g::Union{ZZGenus,Nothing}=nothing; even=true)
  # if given a genus g, this function returns a (list of one) lattice in the genus g that contains L (not necessarily primitively);
  # otherwise, it returns the list of all overlattices M of L up to isometries of M preserving L
  d = ZZ(det(L))
  D = discriminant_group(L)
  idD = hom(D, D, gens(D))
  G,iG = image_in_Oq(L)
  orders = [i for i in divisors(d) if divides(d, i^2)[1]]
  result = ZZLat[]
  for ord in orders 
    #@show ord, D
    b, l, p = is_prime_power_with_data(ord)
    if b && is_elementary(D, p)
      sg = first.(first.(_isotropic_subspaces_representatives_and_stabilizers_elementary(D, iG, valuation(ord, p); do_stab=false)))
    else 
      # slooow
      sg = domain.(first.(_subgroups_orbit_representatives_and_stabilizers(idD, G, ord)))
    end
    for S in sg 
      M = cover(S)
      if !is_integral(M) || (even && !is_even(M))
        continue
      end
      if (g !== nothing)
        em2 = primitive_embeddings(g, M; classification=:first)
        if em2[1] == true
          return [em2[2][1][1]] #returns list of one element for type stability
        end
      else
        push!(result,M)
      end
    end
  end
  return result
end

function root_overlattices(n::Int)
  result = ZZLat[]
  for R in root_lattices(n)
    for S in _overlattice_orbits(R)
      # only add the ones not adding new roots 
      RS = root_sublattice(S)
      if RS == R # this is terribly inefficient
        push!(result,S)
      end 
      if S != RS
        #@show "new"
      end
    end 
  end 
  return result
end


function invariant_function_graph_hash(L::ZZLat; max_size = 6000)
  @assert is_integral(L)
  if rank(L)==0 
    return BigInt(0)
  end
  lb = minimum(L)
  ub = lb + 2*scale(L)
  G = ZZ.(gram_matrix(L))
  # we just care about stuff modulo 2
  n = 0
  kk = GF(2)
  Gk = kk.(G)
  r = rank(L)
  sv = FqMatrix[]
  success = true
  gamma = graph(Undirected, 0)
  _v = zero_matrix(kk, 1, r)
  tmp = zero_matrix(kk, r, 1)
  tmp2 = zero_matrix(kk, 1, 1)
  cv = Hecke._characteristic_vectors(L)
  if length(cv) < max_size
    for v in cv
      n = n+1
      for i in 1:r
        _v[1,i] = v[i]
        tmp[i,1] = v[i]
      end
      Gv = mul!(tmp, Gk, tmp)
      add_vertex!(gamma)
      for (i,w) in enumerate(sv)
        if !iszero(mul!(tmp2, w, Gv))
          add_edge!(gamma, i, n)
        end
      end
      push!(sv, deepcopy(_v))
    end
    return BigInt(canonical_hash(gamma))
  end
  return BigInt(automorphism_group_order(L))
end

function invariant_function_2_4(L::ZZLat; max_size = 10000)  
  v2 = short_vectors(L, 2, Int)
  v4 = short_vectors(L, 4, 4, Int)
  G = Hecke._int_matrix_with_overflow(ZZ.(gram_matrix(L)), ZZ())
  m = MSet{MSet{Int}}()
  #Gi = zero_matrix(ZZ,rank(L),1)
  #v = zero_matrix(ZZ,rank(L),1)
  Gv = zeros(Int, rank(L))
  for (v,_) in v2
    Gv = G*v
    ii = MSet{Int}()
    for (j,_) in v4
      # put the absolute value
      # beause only v or -v is returned by the short vector functions 
      d = abs(dot(j, Gv))
      push!(ii, d)
    end
    push!(m, ii)
  end
  return m
end

function _default_invariant_function(L::ZZLat)
  kn = kissing_number(L)::Int
  rlr, _ = root_lattice_recognition(L)
  R = BigInt
  m = R(ZZ(minimum(L)))
  d = R(ZZ(det(L)))
  ago = R(automorphism_group_order(L))
  if rank(L)>sum(i[2] for i in rlr;init=0)
    t = invariant_function_2_4(L; max_size = 10000)
  else
    t = multiset(multiset([0]))
  end 
  return (m, rlr, kn, ago, d, t)
end


function oscar_invariant_function(L::ZZLat)
  _invariants = Any[]
  _L = L
  while rank(_L) > 12
    M, P, _ = Hecke._shortest_vectors_sublattice(_L; check=false)
    i = index(P,M)
    push!(_invariants, (_default_invariant_function(rescale(P, 1//scale(P); cached=false)),i))
    _L = orthogonal_submodule(_L, P)
  end
  push!(_invariants, invariant_function_graph_hash(_L))
  push!(_invariants, invariant_function_2_4(L; max_size = 10000))
  push!(_invariants, BigInt(automorphism_group_order(L)))
  return Tuple(_invariants)
end


_get_canonical_form(A::ZZMatrix, char_vectors_set::Vector{Matrix{Int}}, canonical_ordering::Vector{Int}) = _get_canonical_form(A, [matrix(ZZ, v) for v in char_vectors_set], canonical_ordering)

function _get_canonical_form(A::ZZMatrix, char_vectors_set::Vector{ZZMatrix}, canonical_ordering::Vector{Int})
  p = length(char_vectors_set)
  filter!(e->e!=p+1 && e!=p+2, canonical_ordering)
  can_char_vectors_set = transpose(matrix(ZZ, reduce(vcat, char_vectors_set[canonical_ordering])))
  _, U = hnf_with_transform(can_char_vectors_set) 
  U_inv = inv(U)
  return transpose(U_inv)*A*U_inv
end

_get_edge_labeled_graph(cv_set::Vector{Matrix{Int}}, gram::ZZMatrix) = _get_edge_labeled_graph(cv_set, Hecke._int_matrix_with_overflow(gram, ZZ(0)))

function _get_edge_labeled_graph(cv_set::Vector{Matrix{Int}}, gram::Matrix{Int})
  p = length(cv_set)
  res_graph = graph(Undirected, p+2)
  max_w = 0
  label!(res_graph, Dict{Tuple{Int, Int}, Int,}(), nothing; name=:edge)
  v_i = Matrix{Int}(undef, 1, number_of_columns(gram))
  t_i = Matrix{Int}(undef, number_of_rows(gram), 1)
  w_i = Matrix{Int}(undef, 1, 1)
  for i = 1:p 
    v_i = AbstractAlgebra.LinearAlgebra.mul!(v_i, cv_set[i], gram)
    for j = i+1:p
      w_i = AbstractAlgebra.LinearAlgebra.mul!(w_i, v_i, AbstractAlgebra.LinearAlgebra.transpose!(t_i, cv_set[j]))
      w = w_i[1]
      max_w = max(w, max_w)
      add_edge!(res_graph, i, j)
      res_graph.edge[i, j] = w
    end
    add_edge!(res_graph, i, p+1)
    w_i = AbstractAlgebra.LinearAlgebra.mul!(w_i, v_i, AbstractAlgebra.LinearAlgebra.transpose!(t_i, cv_set[i]))
    w = w_i[1]
    res_graph.edge[i, p+1] = w
  end
  a = 1+max_w 
  b = a + 1
  for i = 1:p 
    add_edge!(res_graph, i, p+2)
    res_graph.edge[i, p+2] = a
  end
  add_edge!(res_graph, p+1, p+2)
  res_graph.edge[p+1, p+2] = b
  return res_graph
end

"""
    canonical_form(L::ZZLat) -> ZZMatrix
Return the canonical form of ``L``. The form is canonical in the sense, that two isomorphic latticies would have the same canonical form.

We follow ideas of Sikirić, Haensch, Voight and van Woerden [SHVW20](@cite).

!!! note
    We do not give any guarantees that the canonical form stays the same 
    between different versions of Oscar.
"""
function canonical_form(L::ZZLat)
  gram = matrix(ZZ, gram_matrix(L))
  char_vectors_set = Hecke._characteristic_vectors(L)
  # `_characteristic_vectors` returns the characteristic vectors only up to
  # sign; reducing to the fundamental chamber does not commute with the choice
  # of the signs, so we have to take all of them
  char_vectors_set = unique!(append!(char_vectors_set, ZZMatrix[-v for v in char_vectors_set]))
  char_vectors_set = Hecke._reduce_characteristic_vectors(char_vectors_set, L)
  graph = _get_edge_labeled_graph(char_vectors_set, gram) # transform from adjenctcy matrix A to edge-vertex weighted graph Ga, then to edge weighted graph T1(Ga)
  can_order = _canonical_perm(graph; label=:edge) #_canonical_perm uses _edge_label_to_vertex_label themselfs
  return _get_canonical_form(gram, char_vectors_set, can_order)
end
