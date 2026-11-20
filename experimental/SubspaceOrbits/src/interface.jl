################################################################################
#
#  Orbits of matrix groups on subspaces
#
################################################################################

# The data for computing the orbits of G on the k-subspaces of F_p^d.
#
# For 2k > d the computation runs on the dual side: U -> U^perp is a
# bijection between the orbits of G on k-subspaces and the orbits of the
# inverse transposed matrices on (d-k)-subspaces, and the cost of an image
# grows with the number of rows.  Since g -> transpose(inv(g)) is a group
# isomorphism, words in the acting matrices are words in the given ones.
struct SubspaceAction{GT<:MatGroup,K<:LabelAction}
  G::GT
  k::Int
  dual::Bool
  action::K
  # polycyclic generating sequence of a solvable normal subgroup N of G;
  # group elements are kept as GAP matrices, in which stabilizers are built
  pcgs::Vector{GapObj}
  relorders::Vector{Int}
  order_N::ZZRingElem
  # together with N these generate G; empty if GAP is not used, where the
  # action is that of the given generators
  extra::Vector{GapObj}
end

entry_to_int(x::FinFieldElem) = Int(_ffe_to_int(x))
entry_to_int(x::AbstractAlgebra.GFElem) = Int(lift(x))

function int_matrix(M::MatElem)
  return Int[entry_to_int(M[i, j]) for i in 1:nrows(M), j in 1:ncols(M)]
end

# OSCAR has no isomorphism to GAP for the prime fields of AbstractAlgebra,
# so nothing can be asked of GAP about matrix groups over them.
gap_knows(F::Field) = !(F isa AbstractAlgebra.GFField)

# The prime p for a group over F_p.
function field_prime(G::MatGroup)
  q = order(base_ring(G))
  @req is_prime(q) "the group must be defined over a prime field"
  @req q < PACKED_MAX_PRIME "the field must have fewer than $PACKED_MAX_PRIME elements"
  return Int(q)
end

# GAP computes with G via its action on vectors, which is feasible only up
# to about this number of vectors.
const MAX_VECTORS_FOR_GROUP_THEORY = 2^24

# The time GAP needs for the order and the solvable radical of G does not
# depend on the number of subspaces, and can be far longer than applying
# every generator to every subspace.  For lengths alone the radical is
# therefore determined only if it can save at least this many images.
const MIN_IMAGES_SAVED_BY_RADICAL = 2^27

# With `stabilizers = false` the group G is used through its generators
# only, unless its solvable radical can save `min_images_saved` images.
function SubspaceAction(G::MatGroup, k::Int; stabilizers::Bool,
                        min_images_saved::Int=MIN_IMAGES_SAVED_BY_RADICAL)
  p = field_prime(G)
  d = degree(G)
  dual = 2 * k > d
  C = LabelCodec(p, d, dual ? d - k : k)

  pcgs = GapObj[]
  relorders = Int[]
  order_N = ZZ(1)
  extra = GapObj[]
  given = [g for g in unique(gens(G)) if !is_one(g)]
  acting = [matrix(g) for g in given]
  few_vectors = ZZ(p)^d <= MAX_VECTORS_FOR_GROUP_THEORY
  with_gap = few_vectors && gap_knows(base_ring(G))
  if stabilizers
    # they need the order of G, and membership tests in G
    @req few_vectors "stabilizers need at most $MAX_VECTORS_FOR_GROUP_THEORY vectors in F^d"
    @req with_gap "stabilizers need a field as created by `GF`"
  end
  # with the radical about one image per subspace, without it one per generator
  images_saved = (length(given) - 1) * C.D
  if with_gap && (stabilizers || images_saved >= min_images_saved)
    extra = GapObj[GapObj(g) for g in given]
    N = is_solvable(G) ? G : solvable_radical(G)[1]
    gap_pcgs = GAPWrap.Pcgs(GapObj(N))
    pcgs = GapObj[gap_pcgs[i] for i in 1:length(gap_pcgs)]
    relorders = Vector{Int}(GAPWrap.RelativeOrders(gap_pcgs))
    order_N = order(N)
    if order_N == order(G)
      empty!(extra)
    elseif length(extra) > 2
      # the work is proportional to the number of generators
      extra = GapObj[GapObj(g) for g in small_generating_set(G)]
    end
    acting = [matrix(MatGroupElem(G, g)) for g in vcat(pcgs, extra)]
  end

  if dual
    acting = [transpose(inv(g)) for g in acting]
  end
  action = LabelAction(C, Matrix{Int}[int_matrix(g) for g in acting])
  return SubspaceAction(G, k, dual, action, pcgs, relorders, order_N, extra)
end

# The subspace with the given label.
function subspace(A::SubspaceAction, V, lab::Int)
  C = A.action.codec
  M = label_to_matrix(C, lab)
  if A.dual
    M = orthogonal_complement(M, C.p)
  end
  F = base_ring(V)
  return sub(V, [V([F(M[i, j]) for j in 1:ncols(M)]) for i in 1:nrows(M)])[1]
end

################################################################################
#
#  Stabilizers
#
################################################################################

# The element u of N mapping the first point of the orbit to the point at
# 0-based position q of its U_{i+1}-orbit, see `orbits.jl`.
function transversal_element(A::SubspaceAction, state::LabelOrbit, i::Int, q::Int)
  n = length(A.pcgs)
  exps = zeros(Int, n)
  len = prod(view(state.factors, (i + 1):n); init=1)
  for t in (i + 1):n
    len = div(len, state.factors[t])
    exps[t], q = divrem(q, len)
  end
  u = GAPWrap.One(GapObj(A.G))::GapObj
  for t in n:-1:(i + 1)
    exps[t] == 0 && continue
    u = u * A.pcgs[t]^exps[t]
  end
  return u
end

# Generators of the stabilizer in N of the first point of the orbit.  If
# the blocks under U_{i+1} close up after f < (relative order of g_i) steps,
# then pt * g_i^f = pt * u for some u in U_{i+1}, and g_i^f * u^-1 together
# with the stabilizer in U_{i+1} generates the stabilizer in U_i.
function stabilizer_generators_pc(A::SubspaceAction, state::LabelOrbit)
  res = GapObj[]
  for i in 1:length(A.pcgs)
    q = state.closures[i]
    q < 0 && continue
    u = transversal_element(A, state, i, q)
    push!(res, A.pcgs[i]^state.factors[i] * GAPWrap.Inverse(u))
  end
  return res
end

# The product of the h_j on the path from block 0 to block c.
function block_element(A::SubspaceAction, state::LabelOrbit, c::Int)
  path = Int[]
  while c > 0
    c, j = state.parents[c]
    push!(path, j)
  end
  t = GAPWrap.One(GapObj(A.G))::GapObj
  for j in Iterators.reverse(path)
    t = t * A.extra[j]
  end
  return t
end

# Up to this number of edges all of them are enumerated in random order;
# above it edges are drawn at random with repetition.
const MAX_SHUFFLED_EDGES = 2^20

# Extend `gens`, generators of the stabilizer in N of the first point, to
# generators of its stabilizer in G, which has order `target`.
#
# An edge (b, j) of the search over blocks maps the first point pt of the
# orbit, via the first point of block b, to position q of some block c; if
# that is not how block c was found, then t_b * h_j * t_c^-1 * u_q^-1 is a
# Schreier generator of the stabilizer.  Schreier generators of random
# edges are added until the order is reached.
#
# The group they generate is built in the image of the nice monomorphism
# of G, a permutation group, where adding a generator extends the
# stabilizer chain instead of computing a new one.
function extend_stabilizer!(gens::Vector{GapObj}, A::SubspaceAction,
                            state::LabelOrbit, target::ZZRingElem,
                            max_shuffled_edges::Int)
  G = A.G
  K = A.action
  npc = length(A.pcgs)
  m = length(A.extra)
  len0 = state.block_length
  nblocks = div(state.length, len0)
  orbit = view(state.orbit, 1:state.npoints)
  nedges = nblocks * m
  is_tree_edge = falses(nedges)
  for (b, j) in state.parents
    is_tree_edge[b * m + j] = true
  end

  rng = Random.Xoshiro(state.rep)
  shuffled = nedges <= max_shuffled_edges ? Random.randperm(rng, nedges) : Int[]
  next = 0
  nice = GAPWrap.NiceMonomorphism(GapObj(G))
  images = GapObj([GAPWrap.Image(nice, g) for g in gens])
  H = GAPWrap.SubgroupNC(GAPWrap.Image(nice), images)
  while ZZ(GAPWrap.Size(H)) != target
    next += 1
    e = isempty(shuffled) ? rand(rng, 1:nedges) : shuffled[next]
    is_tree_edge[e] && continue
    b, j = divrem(e - 1, m)
    j += 1

    y = image(K, npc + j, Int(orbit[b * len0 + 1])) % eltype(orbit)
    c, q = divrem(findfirst(==(y), orbit)::Int - 1, len0)
    s = block_element(A, state, b) * A.extra[j] *
        GAPWrap.Inverse(transversal_element(A, state, 0, q) * block_element(A, state, c))
    t = GAPWrap.Image(nice, s)
    t in H && continue
    push!(gens, s)
    H = GAPWrap.ClosureGroup(H, t)
    @vprintln :SubspaceOrbits 2 "  $(length(gens)) generators, index $(divexact(target, ZZ(GAPWrap.Size(H))))"
  end
  return gens
end

function point_stabilizer(A::SubspaceAction, state::LabelOrbit;
                          max_shuffled_edges::Int=MAX_SHUFFLED_EDGES)
  G = A.G
  state.length == 1 && return G

  gens = stabilizer_generators_pc(A, state)
  target = divexact(order(G), state.length)
  @vprintln :SubspaceOrbits 1 "orbit of label $(state.rep): length $(state.length), $(div(state.length, state.block_length)) blocks"
  if target != divexact(A.order_N, state.block_length)
    extend_stabilizer!(gens, A, state, target, max_shuffled_edges)
  end
  H = sub(G, elem_type(G)[MatGroupElem(G, g) for g in gens]; check=false)[1]
  set_order(H, target)
  return H
end

################################################################################
#
#  User interface
#
################################################################################

@doc raw"""
    subspace_orbit_representatives_and_lengths(G::MatGroup{<:FinFieldElem}, k::Int)

Return a vector of pairs `(U, n)` where `U` runs over representatives
of the orbits of `G` on the `k`-dimensional subspaces of `F^d`,
and `n` is the length of the orbit of `U`.
Here `G` is a subgroup of `general_linear_group(d, F)`, acting on row vectors
from the right.
The field `F` must be a prime field with fewer than $2^{16}$ elements,
and `F^d` must have at most $2^{36}$ subspaces of dimension `k`.

Each subspace is represented by a single integer during the computation,
see [OBr90](@cite).
The memory needed is one bit per `k`-dimensional subspace of `F^d`,
for the set of subspaces already found,
plus 4 or 8 bytes per subspace in the orbit being computed.

The group `G` is used through its generators only, unless it has many
generators and many subspaces to act on: then its order and its solvable
radical are computed, provided that `F^d` has at most $2^{24}$ vectors.

# Examples
```jldoctest
julia> G = orthogonal_group(1, 4, GF(3))
GO+(4,3)

julia> res = subspace_orbit_representatives_and_lengths(G, 1);

julia> print(sort([n for (U, n) in res]))
[12, 12, 16]
```
"""
function subspace_orbit_representatives_and_lengths(G::MatGroup{<:FinFieldElem}, k::Int)
  d = degree(G)
  @req 0 <= k <= d "k must be between 0 and the degree of the group"
  field_prime(G)
  V = vector_space(base_ring(G), d)
  k == 0 && return [(sub(V, elem_type(V)[])[1], 1)]
  k == d && return [(sub(V, gens(V))[1], 1)]

  A = SubspaceAction(G, k; stabilizers=false)
  labels = Tuple{Int,Int}[]
  label_orbits(A.action, A.relorders; track=false) do state
    push!(labels, (state.rep, state.length))
  end
  return [(subspace(A, V, lab), n) for (lab, n) in labels]
end

@doc raw"""
    subspace_orbit_representatives_and_stabilizers(G::MatGroup{<:FinFieldElem}, k::Int)

Return a vector of pairs `(U, S)` where `U` runs over representatives
of the orbits of `G` on the `k`-dimensional subspaces of `F^d`,
and `S` is the stabilizer of `U` in `G`.
Here `G` is a subgroup of `general_linear_group(d, F)`, acting on row vectors
from the right.
The field `F` must be a prime field with fewer than $2^{16}$ elements,
and `F^d` must have at most $2^{36}$ subspaces of dimension `k`
and, unless `k` is `0` or `d`, at most $2^{24}$ vectors.
The field must be one that GAP can be asked about, as created by `GF`.

The result is of the same form as that of
`orbit_representatives_and_stabilizers`, with possibly different
representatives.
See [`subspace_orbit_representatives_and_lengths`](@ref) for the method.
Unlike that function, this one always computes the order and the
solvable radical of `G`, which can take far longer than the orbits.

# Examples
```jldoctest
julia> G = orthogonal_group(1, 4, GF(3))
GO+(4,3)

julia> res = subspace_orbit_representatives_and_stabilizers(G, 1);

julia> print(sort([index(G, S) for (U, S) in res]))
ZZRingElem[12, 12, 16]
```
"""
function subspace_orbit_representatives_and_stabilizers(G::MatGroup{<:FinFieldElem}, k::Int)
  d = degree(G)
  @req 0 <= k <= d "k must be between 0 and the degree of the group"
  field_prime(G)
  V = vector_space(base_ring(G), d)
  k == 0 && return [(sub(V, elem_type(V)[])[1], G)]
  k == d && return [(sub(V, gens(V))[1], G)]

  A = SubspaceAction(G, k; stabilizers=true)
  res = Tuple{Int,typeof(G)}[]
  label_orbits(A.action, A.relorders; track=true) do state
    push!(res, (state.rep, point_stabilizer(A, state)))
  end
  return [(subspace(A, V, lab), S) for (lab, S) in res]
end
