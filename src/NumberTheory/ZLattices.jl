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
  can_char_vectors_set = transpose(matrix(ZZ, reduce(vcat, char_vectors_set[canonical_ordering])))
  _, U = hnf_with_transform(can_char_vectors_set) 
  U_inv = inv(U)
  return transpose(U_inv)*A*U_inv
end

_canonical_ordering(cv_set::Vector{Matrix{Int}}, gram::ZZMatrix) = _canonical_ordering(cv_set, Hecke._int_matrix_with_overflow(gram, ZZ(0)))

# Return an ordering of the characteristic vectors `cv_set` which depends only
# on the isometry class of the lattice with gram matrix `gram`.
#
# The characteristic vectors together with their inner products form a graph
# whose vertices are colored by the norms and whose edges are labeled by the
# inner products. It is enough to take the pairs meeting a subset which spans
# the rational span of the lattice: a characteristic vector is determined by
# its inner products with such a subset, so every automorphism of the graph is
# still induced by an isometry. We use the classes of equal norm, which are
# invariant, smallest first, until they span. Following section 14 of the nauty
# manual the edge labels are turned into vertex colors: the label of an edge is
# written in binary and the edge is inserted into those of the `nlayers` copies
# of the vertex set which correspond to the bits set. The labels are numbered
# by decreasing frequency, so the most frequent inner product is encoded by
# zero, that is, by no edge at all; this keeps the resulting graph small.
function _canonical_ordering(cv_set::Vector{Matrix{Int}}, gram::Matrix{Int})
  p = length(cv_set)
  cv = reduce(vcat, cv_set)
  A = cv*gram
  norms = [sum(A[i, j]*cv[i, j] for j in 1:size(cv, 2)) for i in 1:p]

  classes = Dict{Int, Vector{Int}}()
  for i in 1:p
    push!(get!(Vector{Int}, classes, norms[i]), i)
  end
  classes = sort!(collect(classes); by = c -> (length(c[2]), c[1]))
  ref = Int[]
  for (_, idx) in classes
    append!(ref, idx)
    rank(matrix(ZZ, cv[ref, :])) == size(cv, 2) && break
  end
  q = length(ref)
  # the characteristic vectors are reordered so that the reference vectors come
  # first; `sigma` translates back
  sigma = append!(copy(ref), setdiff(1:p, ref))
  # rebinding `cv` and `A` here would box them inside `visit_pairs` below
  cvs = cv[sigma, :]
  As = A[sigma, :]
  cls = zeros(Int, p)
  for (c, (_, idx)) in enumerate(classes), i in idx
    cls[i] = c
  end
  col = cls[sigma]

  bs = max(div(2^20, p), 1) # the number of rows treated at once
  # call `f(i, j, w)` for all `i < j` with `i` a reference vector, where `w` is
  # the inner product of the `i`-th and the `j`-th characteristic vector. The
  # inner products are computed in blocks of `bs` rows, so that the whole `q`
  # times `p` matrix of them is never stored.
  function visit_pairs(f::F) where F
    for k in 1:bs:q
      l = min(k+bs-1, q)
      B = As[k:l, :]*transpose(cvs)
      for j in 1:p, i in k:min(l, j-1)
        f(i, j, B[i-k+1, j])
      end
    end
  end

  counts = Dict{Int, Int}()
  visit_pairs((i, j, w) -> counts[w] = get(counts, w, 0) + 1)
  labels = sort!(collect(keys(counts)); by = w -> (-counts[w], w))
  nlayers = length(digits(length(labels)-1; base=2))
  codes = sort(0:2^nlayers-1; by = c -> (count_ones(c), c))
  code = Dict{Int, Int}(w => codes[k] for (k, w) in enumerate(labels))

  colors = [(l-1)*length(classes) + col[i] for l in 1:nlayers for i in 1:p]

  g = graph(Undirected, nlayers*p)
  visit_pairs(function(i, j, w)
    c = code[w]
    while !iszero(c)
      l = trailing_zeros(c)
      add_edge!(g, i + l*p, j + l*p; check=false)
      c &= c-1
    end
  end)
  for l in 1:nlayers-1, i in 1:p
    add_edge!(g, i + (l-1)*p, i + l*p; check=false)
  end
  perm = Polymake._canonical_perm(pm_object(g), Polymake.Array{Int}(colors))
  # the vertices of the first layer come first in the coloring, hence they are
  # exactly the first `p` entries of the canonical labeling
  order = Polymake.to_one_based_indexing(perm)[1:p]
  @assert isperm(order)
  return sigma[order]
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
  L = lll(L) # leaves canonical form unchanged and helps if the basis is badly conditioned
  gram = matrix(ZZ, gram_matrix(L))
  char_vectors_set = Hecke._reduced_characteristic_vectors(L)
  can_order = _canonical_ordering(char_vectors_set, gram)
  return _get_canonical_form(gram, char_vectors_set, can_order)
end
