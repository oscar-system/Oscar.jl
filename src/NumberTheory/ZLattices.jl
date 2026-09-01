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
  cv = characteristic_vectors(L)
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


# return all characteristic vectors up to sign
# unfortunately still to many for a fast graph hash 
# at least in higher rank
# the idea follows https://arxiv.org/pdf/2004.14022
"""
    characteristic_vectors(L::ZZLat) -> Vector{ZZMatrix}
    
Return a set of characteristic vectors of ``L`` up to sign.

We follow ideas of Sikirić, Haensch, Voight and van Woerden [SHVW20](@cite).

!!! note
    We do not give any guarantees that the characteristic vector set stays the same 
    between different versions of Oscar.
"""
function characteristic_vectors(L::ZZLat)
  L = lattice(rational_span(L))
  S1,P1, v1  = Hecke._shortest_vectors_sublattice(L; check=false)
  cvL = v1
  B = coordinates(basis_matrix(S1), P1)
  A = abelian_group(ZZ.(B))
  BS1 = ZZ.(basis_matrix(S1))
  done = []
  for a in A
    -a in done && continue
    iszero(a) && continue
    push!(done,a)
    v = coordinates(a.coeff*basis_matrix(P1), S1)[1,:]
    tmp = [matrix(ZZ, 1, degree(S1), (v -  j)*basis_matrix(S1)) for j in Hecke._closest_vectors(S1, v)[2]]
    append!(cvL, tmp)
  end
  if rank(S1) == rank(L)
    @hassert :Lattice 1 isone(hnf(reduce(vcat, cvL))[1:rank(L),:])
    return cvL
  end
  proj2 = orthogonal_projection(ambient_space(L), basis_matrix(P1))
  # a reduced basis speeds up the closest vector problems in the recursion
  L2 = lll(proj2(L))
  proj1 = orthogonal_projection(ambient_space(L), basis_matrix(L2))
  P_Z = ZZ.(solve(basis_matrix(L2), proj2.matrix;side=:left))    
  ctx = solve_init(P_Z)
  # the differences `w - j` only depend on `w` modulo `P1`, so we solve each
  # closest vector problem only once
  closest = Dict{Vector{QQFieldElem},Vector{Vector{QQFieldElem}}}()
  # recurse 
  for a in characteristic_vectors(L2)
    aL = a*basis_matrix(L2)
    if a*basis_matrix(L2) in L
      @assert rank(L) == ncols(aL)
      push!(cvL, ZZ.(aL))
      continue
    end
    # a vector in L projecting to a
    vL = solve(ctx, a; side=:left)
    w_amb = vL * proj1.matrix
    w_amb == w_amb * proj1.matrix
    w = coordinates(w_amb[1,:], P1)
    if all(isone, denominator.(w))
      push!(cvL, w*basis_matrix(P1))
      continue 
    end
    cv = get!(() -> [w - j for j in Hecke._closest_vectors(P1, w)[2]], closest, [x - floor(x) for x in w])
    tmp = [ZZ.(aL+matrix(QQ, 1, length(w), j) * basis_matrix(P1)) for j in cv]
    append!(cvL, tmp)
  end
  @assert all(rank(L) == ncols(i) for i in cvL)
  @hassert :Lattice 1 isone(hnf(reduce(vcat, cvL))[1:rank(L),:])
  return cvL
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

_reduce_characteristic_vectors(cv_set::Vector{ZZMatrix}, L::ZZLat) = _reduce_characteristic_vectors(_convert_cv_set_to_int(cv_set, matrix(ZZ, gram_matrix(L))), L)
function _reduce_characteristic_vectors(cv_set::Vector{Matrix{Int}}, L::ZZLat)
  R, _, _ = root_lattice_recognition_fundamental(L)
  A = basis_matrix(R)
  gram = matrix(ZZ, gram_matrix(L))
  B_lat = basis_matrix(L)
  A_lat = change_base_ring(ZZ, solve(B_lat, A))
  v_i = Matrix{Int}(undef, 1, number_of_columns(gram))
  t_i = Matrix{Int}(undef, number_of_rows(gram), 1)
  w_i = Matrix{Int}(undef, 1, 1)
  tmp = ZZ(0)
  gram_int = Hecke._int_matrix_with_overflow(gram, tmp)
  A_lat_int = Hecke._int_matrix_with_overflow(A_lat, tmp)
  fundamental_roots = [reshape(A_lat_int[i, :], :, 1) for i in 1:number_of_rows(A_lat_int)]
  res::Vector{Matrix{Int}} = []
  for v in cv_set
    AbstractAlgebra.LinearAlgebra.mul!(v_i, v, gram_int)
    AbstractAlgebra.LinearAlgebra.mul!(w_i, v_i, AbstractAlgebra.LinearAlgebra.transpose!(t_i, v))
    if w_i[1] == 1 || w_i[1] == 2
      continue
    end
    in_chamber = true
    for f_root in fundamental_roots
      AbstractAlgebra.LinearAlgebra.mul!(w_i, v_i, f_root)
      if w_i[1] < 0
        in_chamber = false
        break
      end
    end
    if in_chamber
      push!(res, v)
    end
  end
  for f_root in fundamental_roots
    push!(res, f_root')
  end
  return res
end

function _convert_cv_set_to_int(cv_set::Vector{ZZMatrix}, gram::ZZMatrix)
  tmp = ZZ(0)
  n = ZZ(0)
  tmp1 = zero_matrix(ZZ, 1, number_of_columns(gram))
  tmp2 = zero_matrix(ZZ, number_of_rows(gram), 1)
  tmp3 = zero_matrix(ZZ, 1, 1)
  for v in cv_set
    tmp1 = mul!(tmp1, v, gram)
    tmp3 = mul!(tmp3, tmp1, transpose!(tmp2, v))
    n = max(n, tmp3[1])
  end
  if n-2 < ZZ(typemax(Int)) # we use Cauchy-Schwarz to check if char vector inner products are small enough to be converted to Int. As we need at least w_max+1 and w_max+2 weights further, we need to lower bound by -2.
    cv_set_int = [Hecke._int_matrix_with_overflow(v, tmp) for v in cv_set]
  else 
    throw(OverflowError("The characteristic vectors have to large inner products to be converted to Int."))
  end
  return cv_set_int
end

################################################################################
#
#  Minuscule vectors of root lattices
#
################################################################################

# Let `R` be a root lattice with a fixed fundamental system of roots
# `a_1,...,a_n` and let `c` be a class of the discriminant group `R^vee/R`.
# The vectors of minimal norm of `c` lying in the closed fundamental Weyl
# chamber are called the minuscule vectors of `c`. Since the Weyl group `W(R)`
# acts trivially on `R^vee/R` and the closed fundamental chamber is a
# fundamental domain for its action, they represent the vectors of minimal
# norm of `c` up to the action of `W(R)`. Solving a closest vector problem in
# `R` for a vector of `R^vee` therefore amounts to a table lookup, as long as
# we are only interested in the solutions up to the Weyl group.
#
# We describe a vector `v` of `R^vee` by its weight coordinates
# `g = (v*a_1,...,v*a_n)` in `ZZ^n`. Two such vectors lie in the same class of
# `R^vee/R` if and only if they have the same `g*adj mod d`, where `d` is the
# order of `R^vee/R` and `adj = d*C^-1` for `C` the Cartan matrix, that is to
# say the gram matrix of `a_1,...,a_n`.
struct _MinusculeTable
  d::Int                                                    # order of `R^vee/R`
  adj::Matrix{Int}                                          # `d*C^-1`
  data::Dict{Vector{Int},Tuple{Vector{Int},Rational{Int}}}  # class -> minuscule vector and its norm
end

const _minuscule_tables = Dict{Matrix{Int},Union{_MinusculeTable,Nothing}}()

_minuscule_class(g::AbstractVector{Int}, adj::Matrix{Int}, d::Int) = Int[mod(sum(g[k]*adj[k, j] for k in 1:length(g)), d) for j in 1:length(g)]

# `d` times the norm of the vector with weight coordinates `y`
_minuscule_norm(y::Vector{Int}, adj::Matrix{Int}) = sum(y[i]*adj[i, j]*y[j] for i in 1:length(y), j in 1:length(y); init=0)

# Append to `res` the weight coordinates of all vectors in the closed
# fundamental chamber of norm at most `bound//d`. As the entries of `adj` are
# non-negative, the norm is non-decreasing in each of the (non-negative)
# weight coordinates, which we use to cut the search tree.
function _dominant_weights!(res::Vector{Vector{Int}}, y::Vector{Int}, adj::Matrix{Int}, bound::Int, i::Int)
  if i > length(y)
    push!(res, copy(y))
    return nothing
  end
  k = 0
  while true
    y[i] = k
    _minuscule_norm(y, adj) > bound && break
    _dominant_weights!(res, y, adj, bound, i + 1)
    k += 1
  end
  y[i] = 0
  return nothing
end

# Return the minuscule vectors of the irreducible root lattice with Cartan
# matrix `cartan`, or `nothing` if our assumptions on `cartan` fail. Every
# class of the discriminant group meets the closed fundamental chamber, so
# enumerating the vectors of small norm in there yields the minimal norm of
# each class together with its minuscule vectors. We insist that every class
# has a single minuscule vector, which holds for root lattices of type ADE.
function _minuscule_table(cartan::Matrix{Int})
  n = size(cartan, 1)
  all(i -> cartan[i, i] == 2, 1:n) || return nothing
  dz = det(matrix(ZZ, cartan))
  (dz > 0 && fits(Int, dz)) || return nothing
  d = Int(dz)
  adjq = d*inv(matrix(QQ, cartan))
  adj = Matrix{Int}(undef, n, n)
  for i in 1:n, j in 1:n
    a = adjq[i, j]
    (isone(denominator(a)) && numerator(a) >= 0 && fits(Int, numerator(a))) || return nothing
    adj[i, j] = Int(numerator(a))
  end
  # needed for the enumeration to be complete and to terminate
  all(i -> adj[i, i] > 0, 1:n) || return nothing
  bound = maximum(adj[i, i] for i in 1:n)
  for _ in 1:20
    weights = Vector{Int}[]
    _dominant_weights!(weights, zeros(Int, n), adj, bound, 1)
    norms = Dict{Vector{Int},Rational{Int}}()
    for w in weights
      c = _minuscule_class(w, adj, d)
      nrm = _minuscule_norm(w, adj)//d
      norms[c] = min(get(norms, c, nrm), nrm)
    end
    if length(norms) == d
      data = Dict{Vector{Int},Tuple{Vector{Int},Rational{Int}}}()
      for w in weights
        c = _minuscule_class(w, adj, d)
        _minuscule_norm(w, adj)//d == norms[c] || continue
        haskey(data, c) && return nothing
        data[c] = (w, norms[c])
      end
      return _MinusculeTable(d, adj, data)
    end
    bound *= 2
  end
  return nothing
end

# The tables of the irreducible components of a root lattice with Cartan
# matrix `cartan`, the components being given by `ranges`
function _minuscule_tables_of(cartan::ZZMatrix, ranges::Vector{UnitRange{Int}})
  res = Tuple{UnitRange{Int},_MinusculeTable}[]
  for r in ranges
    all(x -> -2 <= x <= 2, (cartan[i, j] for i in r for j in r)) || return nothing
    block = Int[Int(cartan[i, j]) for i in r, j in r]
    t = get!(() -> _minuscule_table(block), _minuscule_tables, block)
    t === nothing && return nothing
    push!(res, (r, t))
  end
  return res
end

# The weight coordinates of the minuscule vector of the class of the vector
# with weight coordinates `g`, together with its norm. Both the closed
# fundamental chamber and the discriminant group are direct products over the
# irreducible components, so that we may work component by component.
function _minuscule_vector(D::Vector{Tuple{UnitRange{Int},_MinusculeTable}}, g::ZZMatrix)
  gd = zeros(Int, ncols(g))
  nrm = zero(Rational{Int})
  for (r, t) in D
    w, s = t.data[_minuscule_class(Int[Int(mod(g[1, j], t.d)) for j in r], t.adj, t.d)]
    gd[r] = w
    nrm += s
  end
  return gd, nrm
end

# The vector of minimal norm of the coset `v + R` lying in the closed
# fundamental chamber, where `R` is the root lattice with fundamental roots
# the rows of `roots`, `g` are the weight coordinates of the vector `v` and
# `gd` those of the minuscule vector of its class. The difference lies in `R`
# and has coordinates `(gd - g)*C^-1` with respect to the fundamental roots.
function _minuscule_translate(v::ZZMatrix, g::ZZMatrix, gd::Vector{Int}, D::Vector{Tuple{UnitRange{Int},_MinusculeTable}}, roots::ZZMatrix)
  z = zero_matrix(ZZ, 1, nrows(roots))
  for (r, t) in D, (jj, j) in enumerate(r)
    z[1, j] = divexact(sum((gd[k] - g[1, k])*t.adj[kk, jj] for (kk, k) in enumerate(r)), t.d)
  end
  return v + z*roots
end

################################################################################
#
#  Reduced characteristic vectors
#
################################################################################

# Return the reduced characteristic vector set of `L`, that is to say the
# fundamental roots of `L` together with the characteristic vectors of norm
# different from 1 and 2 lying in the closed fundamental Weyl chamber.
#
# Return `nothing` unless the vectors of minimal norm of `L` are roots. If
# they are, then the sublattice they span is the root sublattice `R` of `L`,
# and all closest vector problems occurring in `characteristic_vectors` are
# closest vector problems in `R` for vectors of `R^vee`. Their solutions are
# the minuscule vectors of `R` up to the action of the Weyl group of `R`, and
# these are all we need here: the Weyl group is contained in the orthogonal
# group of `L`, so keeping only the solutions in the closed fundamental
# chamber loses no information. Note that all the vectors we produce have norm
# bigger than 2: they are non-zero and none of them is a root, since a root
# lies in `R` and has trivial class in `R^vee/R`.
function _reduced_characteristic_vectors_with_roots(L::ZZLat; check=true)
  #check && ((is_positive_definite(L) && is_integral(L) && minimum(L) == 2) || return nothing)
  n = rank(L)
  gram = change_base_ring(ZZ, gram_matrix(L))
  # the fundamental roots of `L`, in the coordinates of `L` and grouped into
  # the irreducible components of the root sublattice
  _, components = Hecke._root_lattice_recognition_fundamental(L)
  roots = reduce(vcat, components; init=zero_matrix(ZZ, 0, n))
  nr = nrows(roots)
  ranges = UnitRange{Int}[]
  k = 0
  for c in components
    push!(ranges, (k+1):(k+nrows(c)))
    k += nrows(c)
  end
  cartan = roots*gram*transpose(roots)
  D = _minuscule_tables_of(cartan, ranges)
  D === nothing && return nothing
  # `v*gram_roots` are the weight coordinates of the vector `v` of `L`
  gram_roots = gram*transpose(roots)
  res = ZZMatrix[roots[i:i, :] for i in 1:nr]

  # the vectors of minimal norm of the non-trivial cosets of `P/R`, where `P`
  # is the primitive closure of `R` in `L`, are characteristic vectors
  BP = saturate(roots)
  cosets = ZZMatrix[zero_matrix(ZZ, 1, n)]
  # the square of the index of `R` in `P`
  if !isone(divexact(det(cartan), det(BP*gram*transpose(BP))))
    for a in abelian_group(solve(BP, roots; side=:left))
      is_zero(a) && continue
      push!(cosets, a.coeff*BP)
    end
  end
  for v in cosets[2:end]
    g = v*gram_roots
    push!(res, _minuscule_translate(v, g, _minuscule_vector(D, g)[1], D, roots))
  end
  nr == n && return _convert_cv_set_to_int(res, gram)

  # the remaining characteristic vectors are the lifts of minimal norm of the
  # characteristic vectors of the projection `L2` of `L` to the orthogonal
  # complement of `R`
  L0 = lattice(rational_span(L))
  P = lattice(ambient_space(L0), change_base_ring(QQ, BP); isbasis=true, check=false)
  proj = orthogonal_projection(ambient_space(L0), basis_matrix(P))
  L2 = proj(L0)
  PZ = change_base_ring(ZZ, solve(basis_matrix(L2), proj.matrix; side=:left))
  cv2 = characteristic_vectors(L2)
  # `characteristic_vectors` returns the characteristic vectors only up to
  # sign; we need all of them for the result to be canonical
  cv2 = unique!(append!(cv2, ZZMatrix[-a for a in cv2]))
  ctx = solve_init(PZ)
  gs = Vector{ZZMatrix}(undef, length(cosets))
  gds = Vector{Vector{Int}}(undef, length(cosets))
  ns = Vector{Rational{Int}}(undef, length(cosets))
  for a in cv2
    # the lifts of `a` are the vectors of `vL + P`; their norm is the norm of
    # `a` plus the norm of the part in `P`, so it is minimal exactly for the
    # cosets of `P/R` with the shortest minuscule vectors
    vL = solve(ctx, a; side=:left)
    for (i, c) in enumerate(cosets)
      gs[i] = (vL + c)*gram_roots
      gds[i], ns[i] = _minuscule_vector(D, gs[i])
    end
    best = minimum(ns)
    for (i, c) in enumerate(cosets)
      ns[i] == best || continue
      push!(res, _minuscule_translate(vL + c, gs[i], gds[i], D, roots))
    end
  end
  return _convert_cv_set_to_int(res, gram)
end

# Return the fundamental roots of `L` together with the characteristic vectors
# of norm different from 1 and 2 lying in the closed fundamental Weyl chamber.
function _reduced_characteristic_vectors(L::ZZLat)
  if gram_matrix(L)[1,1] < 0
    L = rescale(L, -1)
  end
  res = _reduced_characteristic_vectors_with_roots(L)
  res === nothing || return res
  cv = characteristic_vectors(L)
  # `characteristic_vectors` returns the characteristic vectors only up to
  # sign; we need all of them for the result to be canonical
  cv = unique!(append!(cv, ZZMatrix[-v for v in cv]))
  return _reduce_characteristic_vectors(cv, L)
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
  # the computations below are much faster with a reduced basis, and the
  # result does not depend on the chosen basis
  L = lll(L; _is_definite=true)
  gram = matrix(ZZ, gram_matrix(L))
  char_vectors_set = _reduced_characteristic_vectors(L)
  graph = _get_edge_labeled_graph(char_vectors_set, gram) # transform from adjenctcy matrix A to edge-vertex weighted graph Ga, then to edge weighted graph T1(Ga)
  can_order = _canonical_perm(graph; label=:edge) #_canonical_perm uses _edge_label_to_vertex_label themselfs
  return _get_canonical_form(gram, char_vectors_set, can_order)
end
