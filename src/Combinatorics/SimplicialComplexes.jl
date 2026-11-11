import Oscar: Polymake, pm_object


################################################################################
##  Constructing
################################################################################

mutable struct SimplicialComplex
    pm_simplicialcomplex::Polymake.BigObject
    face_cache::Dict{Int,Vector{Set{Int}}}
    SimplicialComplex(bo::Polymake.BigObject) = new(bo)
end

pm_object(K::SimplicialComplex) = K.pm_simplicialcomplex


@doc raw"""
    simplical_complex(generators::Union{Vector{Vector{Int}}, Vector{Set{Int}}})

Construct an abstract simplicial complex from a set of faces.
While arbitrary non-negative integers are allowed as vertices, they will be relabeled to consecutive integers starting at 1.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]])
Abstract simplicial complex of dimension 2 on 4 vertices

julia> G = complete_bipartite_graph(2,3)
Undirected graph with 5 nodes and the following edges:
(3, 1)(3, 2)(4, 1)(4, 2)(5, 1)(5, 2)

julia> K = simplicial_complex(G)
Abstract simplicial complex of dimension 1 on 5 vertices
```

Simplicial complex comprising the empty set only:
```jldoctest
julia> empty = simplicial_complex(Vector{Set{Int}}([]))
Abstract simplicial complex of dimension -1 on 0 vertices
```

The original vertices can be recovered:
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> L = simplicial_complex([[0,2,17],[2,17,90]]);

julia> facets(L)
2-element Vector{Set{Int64}}:
 Set([2, 3, 1])
 Set([4, 2, 3])

julia> vertexindices(L)
4-element Vector{Int64}:
  0
  2
 17
 90
```
"""
function simplicial_complex(generators::Union{AbstractVector{<:AbstractVector{<:Base.Integer}}, AbstractVector{<:AbstractSet{<:Base.Integer}}})
  K = Polymake.topaz.SimplicialComplex(INPUT_FACES=generators)
  SimplicialComplex(K)
end

function simplicial_complex(generators::Union{Set{<:AbstractVector{<:Base.Integer}}, Set{<:AbstractSet{<:Base.Integer}}})
  return simplicial_complex(collect(generators))
end

function simplicial_complex(generators::IncidenceMatrix)
  K = Polymake.@convert_to Array{Set} Polymake.common.rows(generators)
  simplicial_complex(K)
end

function simplicial_complex(G::Graph)
  IM = incidence_matrix(G)
  edges_as_rows = IncidenceMatrix(transpose(IM))
  simplicial_complex(edges_as_rows)
end

# more efficient UNEXPORTED+UNDOCUMENTED version, which requires consecutive vertices, and facets as generators;
# will produce errors / segfaults or worse if used improperly
function _SimplicialComplex(generators::Union{Vector{Vector{Int}}, Vector{Set{Int}}})
  K = Polymake.topaz.SimplicialComplex(FACETS=generators)
  SimplicialComplex(K)
end
    
################################################################################
##  Auxiliary
################################################################################

function vertexindices(K::SimplicialComplex) 
    bigobject = pm_object(K)
    if Polymake.exists(bigobject,"VERTEX_INDICES")
        return Vector{Int}(bigobject.VERTEX_INDICES)
    else
        return Vector{Int}(1:bigobject.N_VERTICES)
    end
end

function _convert_finitely_generated_abelian_group(A::Polymake.Array{Polymake.HomologyGroup{Polymake.Integer}},
                                                   h_index::Int)
  # we return non-reduced homology and cohomology
  B = A[h_index + 1] # index shift
  betti_number = is_zero(h_index) ? Polymake.betti_number(B) + 1 : Polymake.betti_number(B)
  vec = zeros(Int, betti_number)
  torsion_i = Polymake.torsion(B)
  for (p,k) in torsion_i
    append!(vec, fill(p,k))
  end
  return abelian_group(vec)
end

################################################################################
##  Properties
################################################################################

@doc raw"""
    n_vertices(K::SimplicialComplex)

Return the number of vertices of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> n_vertices(torus())
7
```
"""
n_vertices(K::SimplicialComplex) = pm_object(K).N_VERTICES::Int

@doc raw"""
    n_facets(K::SimplicialComplex)

Return the number of facets of the abstract simplicial complex `K`.
"""
n_facets(K::SimplicialComplex) = pm_object(K).N_FACETS::Int

@doc raw"""
    facets(K::SimplicialComplex)

Return the maximal (by inclusion) faces of the abstract simplicial complex `K`.
"""
function facets(K::SimplicialComplex)
    bigobject = pm_object(K)
    the_facets = Polymake.to_one_based_indexing(bigobject.FACETS)
    return Vector{Set{Int}}(the_facets)
end

function _get_face_cache(K::SimplicialComplex)
  if !isdefined(K, :face_cache)
    K.face_cache = Dict{Int,Vector{Set{Int}}}()
  end
  return K.face_cache
end

@doc raw"""
    faces(K::SimplicialComplex [, dim::Int])

Return the faces of the abstract simplicial complex `K`, passing `dim` as the second argument returns faces of dimension `dim` (sets of size `dim + 1`).
"""
function faces(K::SimplicialComplex, d::Int)
  @req -1 <= d <= dim(K) "Dimension must be between -1 and dim(K)"
  fc = _get_face_cache(K)
  return get!(fc, d) do
    po = face_poset(K)
    Set.(data.(elements_of_rank(po, d + 1)))
  end
end

function faces(K::SimplicialComplex)
  return reduce(vcat, faces(K,i) for i in -1:dim(K))
end

@doc raw"""
    dim(K::SimplicialComplex)

Return the dimension of the abstract simplicial complex `K`.
"""
dim(K::SimplicialComplex) = pm_object(K).DIM::Int

@doc raw"""
    f_vector(K::SimplicialComplex)

Return the face vector (number of faces per dimension) of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> f_vector(torus())
3-element Vector{Int64}:
  7
 21
 14
```
"""
f_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).F_VECTOR)

@doc raw"""
    h_vector(K::SimplicialComplex)

Return the h-vector of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> h_vector(torus())
4-element Vector{Int64}:
  1
  4
 10
 -1
```
"""
h_vector(K::SimplicialComplex) = Vector{Int}(pm_object(K).H_VECTOR)

@doc raw"""
    betti_numbers([R::Union{<: Field, ZZRing}=ZZ,] K::SimplicialComplex)

Return the reduced Betti numbers of the abstract simplicial complex `K`.
Defaults to computing Betti numbers over `ZZ`, otherwise computes the Betti numbers over the ring `R`.

# Examples
```jldoctest
julia> betti_numbers(klein_bottle())
3-element Vector{Int64}:
 0
 1
 0

julia> betti_numbers(GF(2), klein_bottle())
3-element Vector{Int64}:
 0
 2
 1
```
"""
betti_numbers(K::SimplicialComplex) = Vector{Int}(Polymake.topaz.betti_numbers(pm_object(K)))

# Boundary matrix of `K` over `R` from `d`-faces (rows) to `(d-1)`-faces (columns), both in the
# order of `faces`, with the conventions of `Polymake.topaz.boundary_matrix`: a column of ones for
# `d <= 0` and a zero row for `d > dim(K)`.
# ponytail: replaces topaz boundary_matrix, which is wrong for non-pure complexes
# (stray entries in row 0 from facets adjacent to the top node of the Hasse diagram);
# switch back once the polymake fix is released.
function _boundary_matrix(R::Ring, K::SimplicialComplex, d::Int)
  d > dim(K) && return zero_matrix(R, 1, length(faces(K, d - 1)))
  F = faces(K, d)
  d <= 0 && return ones_matrix(R, length(F), 1)
  idx = Dict(s => j for (j, s) in enumerate(faces(K, d - 1)))
  M = zero_matrix(R, length(F), length(idx))
  for (i, f) in enumerate(F), (p, v) in enumerate(sort!(collect(f)))
    M[i, idx[setdiff(f, v)]] = isodd(p) ? one(R) : -one(R)
  end
  return M
end

function betti_numbers(R::Union{<:Field, ZZRing}, K::SimplicialComplex)
  c = characteristic(R)
  iszero(c) && return betti_numbers(K)
  b = Int[]
  boundary_m = _boundary_matrix(R, K, 0)
  im_dim = 1
  for k = 1:dim(K) + 1
    ker_dim = size(boundary_m)[1] - im_dim
    boundary_m = _boundary_matrix(R, K, k)
    im_dim = rank(boundary_m)
    
    push!(b, ker_dim - im_dim)
  end
  return b
end

@doc raw"""
    euler_characteristic(K::SimplicialComplex)

Return the reduced Euler characteristic of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> euler_characteristic(complex_projective_plane())
2
```
"""
euler_characteristic(K::SimplicialComplex) = pm_object(K).EULER_CHARACTERISTIC::Int

@doc raw"""
    homology(K::SimplicialComplex, i::Int)

Return `i`-th integral homology group of `K`.

# Examples
```jldoctest
julia> [ homology(real_projective_plane(), i) for i in [0,1,2] ]
3-element Vector{FinGenAbGroup}:
 Z
 Z/2
 Z/1
```
"""
homology(K::SimplicialComplex, i::Int) = _convert_finitely_generated_abelian_group(pm_object(K).HOMOLOGY, i)

@doc raw"""
    cohomology(K::SimplicialComplex, i::Int)

Return `i`-th integral cohomology group of `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[0,1],[1,2],[0,2]]);

julia> cohomology(K,1)
Z
```
"""
cohomology(K::SimplicialComplex, i::Int) = _convert_finitely_generated_abelian_group(pm_object(K).COHOMOLOGY, i)

@doc raw"""
    minimal_nonfaces(K::SimplicialComplex)

Return the minimal non-faces of the abstract simplicial complex `K`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> minimal_nonfaces(K)
1-element Vector{Set{Int64}}:
 Set([4, 1])
```
"""
minimal_nonfaces(K::SimplicialComplex) = minimal_nonfaces(Vector{Set{Int}}, K)
function minimal_nonfaces(::Type{Vector{Set{Int}}}, K::SimplicialComplex)
    I = minimal_nonfaces(IncidenceMatrix, K)
    return Vector{Set{Int}}([Polymake.row(I,i) for i in 1:Polymake.nrows(I)])
end
function minimal_nonfaces(::Type{IncidenceMatrix}, K::SimplicialComplex)
    # the following line must stay to ensure polymake uses the correct algorithm for the non-faces
    nv = n_vertices(K)
    m = pm_object(K).MINIMAL_NON_FACES
    # fix column number (see #1440) until this is fixed in polymake
    if size(m, 2) < nv
      resize!(m, size(m, 1), nv)
    end
    return m
end

@doc raw"""
    alexander_dual(K::SimplicialComplex)

Return the Alexander dual of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> alexander_dual(K)
Abstract simplicial complex of dimension 1 on 2 vertices
```
"""
alexander_dual(K::SimplicialComplex) = SimplicialComplex(Polymake.topaz.alexander_dual(pm_object(K)))

@doc raw"""
    stanley_reisner_ideal(K::SimplicialComplex)

Return the Stanley-Reisner ideal of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> stanley_reisner_ideal(real_projective_plane())
Ideal generated by
  x1*x2*x3
  x1*x2*x4
  x1*x5*x6
  x2*x5*x6
  x1*x3*x6
  x1*x4*x5
  x3*x4*x5
  x3*x4*x6
  x2*x3*x5
  x2*x4*x6
```
"""
function stanley_reisner_ideal(K::SimplicialComplex)
    n = n_vertices(K)
    R, _ = polynomial_ring(QQ, n, cached=false)
    return stanley_reisner_ideal(R, K)
end

@doc raw"""
    stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)

Return the Stanley-Reisner ideal of the abstract simplicial complex `K`, in the given ring `R`.

# Examples
```jldoctest
julia> R, _ = QQ[:a, :b, :c, :d, :e, :f];

julia> stanley_reisner_ideal(R, real_projective_plane())
Ideal generated by
  a*b*c
  a*b*d
  a*e*f
  b*e*f
  a*c*f
  a*d*e
  c*d*e
  c*d*f
  b*c*e
  b*d*f
```
"""
function stanley_reisner_ideal(R::MPolyRing, K::SimplicialComplex)
  mnf = minimal_nonfaces(IncidenceMatrix, K)
  gens = [ R([1], [Vector{Int}(mnf[i,:])]) for i in 1:Polymake.nrows(mnf) ]
  # currently no way to set as a universal groebner basis
  I = IdealGens(R, gens, default_ordering(R); isGB=true, isReduced=true)
  return ideal(I)
end

@doc raw"""
    stanley_reisner_ring(K::SimplicialComplex)

Return the Stanley-Reisner ring of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> stanley_reisner_ring(K)
(Quotient of multivariate polynomial ring by ideal (x1*x4), Map: multivariate polynomial ring -> quotient of multivariate polynomial ring)
```
"""
function stanley_reisner_ring(K::SimplicialComplex)
    n = n_vertices(K)
    R, _ = polynomial_ring(QQ, n, cached=false)
    return stanley_reisner_ring(R, K)
end

@doc raw"""
    stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex)

Return the Stanley-Reisner ring of the abstract simplicial complex `K`, as a quotient of a given ring `R`.

# Examples
```jldoctest
julia>  R, _ = ZZ[:a, :b, :c, :d, :e, :f];

julia> stanley_reisner_ring(R, real_projective_plane())
(Quotient of multivariate polynomial ring by ideal (a*b*c, a*b*d, a*e*f, b*e*f, a*c*f, a*d*e, c*d*e, c*d*f, b*c*e, b*d*f), Map: R -> quotient of multivariate polynomial ring)
```
"""
stanley_reisner_ring(R::MPolyRing, K::SimplicialComplex) = quo(R, stanley_reisner_ideal(R, K))

_support(m) = findall(>(0), exponent_vector(m, 1))

function stanley_reisner_complex(monomials::Vector{<:MPolyRingElem})
  n = ngens(parent(first(monomials)))
  delta_star = simplicial_complex([setdiff(1:n, _support(mon)) for mon in monomials])
  return alexander_dual(delta_star)
end

# Generators of the upper Koszul simplicial complex of the monomial ideal with generators
# of exponent vectors `gens` in degree `b`: the squarefree tau with x^(b - tau) in the ideal.
# Such a tau lies below {j : g[j] < b[j]} for a generator g dividing x^b, and conversely.
# No generators means the void complex, a single empty generator the complex {emptyset}.
_upper_koszul_generators(gens::Vector{Vector{Int}}, b::Vector{Int}) =
  [findall(g .< b) for g in gens if all(g .<= b)]

# squarefree monomial -> bitmask (assumes ≤ 64 variables)
function to_mask(m)
  b = zero(UInt64)
  for (i, e) in enumerate(exponent_vector(m, 1))
    e > 0 && (b |= one(UInt64) << (i - 1))
  end
  b
end

# keep only minimal generators (drop g if some other h divides it)
minimal_masks(gs) = (u = unique(gs); filter(g -> !any(h -> h != g && h & g == h, u), u))

# L as a sorted vector of masks, including 0 (the bottom element, lcm of nothing)
function _lcm_lattice_squarefree(I::Vector{<:MPolyRingElem})
  L = UInt64[0]
  for g in minimal_masks(to_mask.(I))
    L = unique!(sort!(vcat(L, L .| g)))
  end
  L
end

# Mayer-Vietoris tree of the squarefree monomial ideal generated by `I` (as bitmasks).
# Node <g_1, ..., g_k> at dimension d has left child <g_1, ..., g_{k-1}> at d and
# right child <lcm(g_1, g_k), ..., lcm(g_{k-1}, g_k)> at d + 1. Relevant nodes are the root
# and the right children; walking the left spine of a relevant node visits the right child
# of every pivot g_k, so only relevant nodes go on the stack.
# Returns (d, mask) => number of relevant nodes at dimension d having mask as a generator.
# This bounds beta_{d, mask}(I) from above, with equality if mask occurs in only one d.
function _mayer_vietoris_tree(I::Vector{<:MPolyRingElem})
  counts = Dict{Tuple{Int, UInt64}, Int}()
  stack = [(minimal_masks(to_mask.(I)), 0)]
  while !isempty(stack)
    gens, d = pop!(stack)
    for (k, g) in enumerate(gens)
      counts[(d, g)] = get(counts, (d, g), 0) + 1
      k > 1 && push!(stack, (minimal_masks([h | g for h in gens[1:k-1]]), d + 1))
    end
  end
  counts
end

function _hochster_betti_numbers_squarefree(monomials::Vector{<:MPolyRingElem})
  R = parent(first(monomials))
  n = ngens(R)
  F = coefficient_ring(R)
  T = _mayer_vietoris_tree(monomials)
  vec(m) = digits(Int, m; base = 2, pad = n)

  # multidegree -> dimensions it occurs in the tree
  dims = Dict{UInt64, Vector{Int}}()
  for (d, m) in keys(T)
    push!(get!(dims, m, Int[]), d)
  end

  gens = minimal_masks(to_mask.(monomials))
  B = Dict{Tuple{Int, Vector{Int}}, Int}()
  ambiguous = UInt64[]
  for (m, ds) in dims
    if length(ds) == 1   # tree count is exact
      B[(only(ds), vec(m))] = T[(only(ds), m)]
    else
      push!(ambiguous, m)
    end
  end

  # Hochster: beta_{d, m}(I) = dim H~_{|m| - d - 2}(Delta|_m), where Delta|_m consists of the
  # subsets of m containing no generator (Alexander dual to the upper Koszul complex K^m, and
  # much smaller when there are many generators). m <= m' implies Delta|_m <= Delta|_m', so a
  # chain of multidegrees is a filtration and one persistence computation gives the homology
  # of every Delta|_m on it.
  # ponytail: greedy chain cover extending only chains whose top is m minus one variable
  # (a dictionary lookup per variable), Dilworth (bipartite matching) gives the minimum one
  chains = Vector{UInt64}[]
  tops = Dict{UInt64, Int}()   # top of a chain -> its index
  for m in sort!(ambiguous; by = count_ones)
    c, t = 0, m
    while t != 0 && c == 0
      c = pop!(tops, m & ~(one(UInt64) << trailing_zeros(t)), 0)
      t &= t - 1
    end
    c == 0 ? (push!(chains, [m]); c = length(chains)) : push!(chains[c], m)
    tops[m] = c
  end
  @showprogress for chain in chains
    qmax = maximum(count_ones(m) - minimum(dims[m]) - 2 for m in chain)
    betti = _persistent_reduced_betti(F, gens, chain, qmax)
    for (i, m) in enumerate(chain), d in dims[m]
      h = betti(i, count_ones(m) - d - 2)
      iszero(h) || (B[(d, vec(m))] = h)
    end
  end
  return B
end

# Faces of dimension <= p of Delta|_top for top = chain[end] (the subsets containing no
# generator), as bitmasks, with birth the first i such that the face lies in chain[i].
function _restriction_faces(gens::Vector{UInt64}, chain::Vector{UInt64}, p::Int)
  top = chain[end]
  gs = filter(g -> g & top == g, gens)
  birth = Dict{UInt64, Int}()
  stack = [(zero(UInt64), top)]   # (face, vertices above its largest one)
  while !isempty(stack)
    s, cand = pop!(stack)
    birth[s] = findfirst(m -> s & m == s, chain)
    count_ones(s) > p && continue
    while cand != 0
      t = s | (one(UInt64) << trailing_zeros(cand))
      cand &= cand - 1
      any(g -> g & t == g, gs) || push!(stack, (t, cand))
    end
  end
  return birth
end

# Persistent reduced (co)homology over the field F of the filtration
# Delta|_{chain[1]} <= Delta|_{chain[2]} <= ... in dimensions <= qmax.
# Returns (i, q) -> dim H~_q(Delta|_{chain[i]}; F).
function _persistent_reduced_betti(F::Field, gens::Vector{UInt64}, chain::Vector{UInt64}, qmax::Int)
  birth = _restriction_faces(gens, chain, qmax + 1)
  order = sort!(collect(keys(birth)); by = s -> (birth[s], count_ones(s)))
  idx = Dict(s => j for (j, s) in enumerate(order))

  # Persistent cohomology (as in Ripser): reduce the coboundary columns in reverse filtration
  # order, dimensions bottom-up, with clearing: a face which is the pivot of a column one
  # dimension down is negative and its own column is skipped. The persistence pairs are the
  # same as for homology. Row positions are reversed so that the pivot (the earliest cofacet)
  # is the last entry. The empty face is included, so homology is reduced.
  N = length(order)
  top = chain[end]
  pivot = Dict{Int, Int}()   # pivot position -> column having it
  cols = Dict{Int, SRow{elem_type(F)}}()
  death = zeros(Int, N)      # stage at which the class born with this face dies, 0 = never
  negative = falses(N)
  for c in 0:qmax + 1, j in N:-1:1
    s = order[j]
    (count_ones(s) == c && !negative[j]) || continue
    cobd = Tuple{Int, Int}[]
    t = top & ~s
    while t != 0
      v = trailing_zeros(t)
      t &= t - 1
      jt = get(idx, s | (one(UInt64) << v), 0)
      jt == 0 && continue
      # sign of s in the boundary of s + v: (-1)^(number of elements of s below v)
      push!(cobd, (N + 1 - jt, iseven(count_ones(s & ((one(UInt64) << v) - 1))) ? 1 : -1))
    end
    col = sparse_row(F, cobd)
    while !iszero(col)
      k = get(pivot, col.pos[end], 0)
      k == 0 && break
      Hecke.add_scaled_row!(cols[k], col, -col.values[end] // cols[k].values[end])
    end
    iszero(col) && continue
    jt = N + 1 - col.pos[end]
    pivot[col.pos[end]] = j
    cols[j] = col
    death[j] = birth[order[jt]]
    negative[jt] = true
  end

  positive = .!negative
  return (i, q) -> count(j -> positive[j] && count_ones(order[j]) == q + 1 && birth[order[j]] <= i &&
                              (death[j] == 0 || death[j] > i), eachindex(order))
end

function _dual_hochster_betti_numbers_squarefree(monomials::Vector{<:MPolyRingElem})
  R = parent(first(monomials))
  n = ngens(R)
  delta_star = simplicial_complex([setdiff(1:n, _support(mon)) for mon in monomials])

  full = n == 64 ? typemax(UInt64) : (one(UInt64) << n) - 1
  betti_table = Dict{MPolyRingElem, Vector{Int}}()
  @showprogress for sigma in Oscar._lcm_lattice_squarefree(monomials)
    sigma == 0 && continue                                         # bottom element, not a degree
    mon = monomial(R, digits(Int, sigma; base = 2, pad = n))       # key is sigma, not sigma_bar
    if mon in monomials                                            # generator: link is {∅}
      betti_table[mon] = [1]
      continue
    end
    sigma_bar = ~sigma & full
    link = link_subcomplex(delta_star, findall(isone, digits(sigma_bar; base = 2)))  # link of the complement
    betti_table[mon] = betti_numbers(link)
  end
  return betti_table
end

@doc raw"""
    hochster_betti_numbers(monomials::Vector{<:MPolyRingElem})

Return the nonzero multigraded Betti numbers of the ideal $I$ generated by `monomials`,
as a dictionary sending $(i, b)$ to $\beta_{i, b}(I)$, computed by Hochster's formula
```math
\beta_{i, b}(I) = \dim_k \tilde{H}_{i - 1}(K^b(I); k),
```
where $k$ is the coefficient field and
$K^b(I) = \{\text{squarefree } \tau : x^{b - \tau} \in I\}$ is the upper Koszul
simplicial complex. Only degrees $b$ in the lcm lattice of the generators can contribute.
Note that $\beta_{i, b}(I) = \beta_{i + 1, b}(S/I)$.

The monomials do not have to be squarefree or minimal generators. Their coefficients are
ignored.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> B = hochster_betti_numbers([x^2, x*y, y^2]);

julia> sort(collect(B))
5-element Vector{Pair{Tuple{Int64, Vector{Int64}}, Int64}}:
 (0, [0, 2]) => 1
 (0, [1, 1]) => 1
 (0, [2, 0]) => 1
 (1, [1, 2]) => 1
 (1, [2, 1]) => 1
```
"""
function hochster_betti_numbers(monomials::Vector{<:MPolyRingElem})
  @req !isempty(monomials) "The list of monomials must not be empty"
  @req all(is_monomial, monomials) "All elements must be monomials"
  F = coefficient_ring(parent(monomials[1]))
  @req F isa Field "The coefficient ring must be a field"
  gens = unique!([Int.(leading_exponent_vector(m)) for m in monomials])

  lcm_lattice = Set(gens)
  new_degrees = copy(gens)
  while !isempty(new_degrees)
    new_degrees = [d for d in unique!([max.(a, g) for a in new_degrees for g in gens]) if !(d in lcm_lattice)]
    union!(lcm_lattice, new_degrees)
  end

  B = Dict{Tuple{Int, Vector{Int}}, Int}()
  for b in lcm_lattice
    faces = filter(!isempty, _upper_koszul_generators(gens, b))
    # b is in the lcm lattice, so some generator divides x^b: the complex is not void
    if isempty(faces)
      B[(0, b)] = 1
      continue
    end
    for (i, d) in enumerate(betti_numbers(F, simplicial_complex(faces)))
      iszero(d) || (B[(i, b)] = d)
    end
  end
  return B
end

################################################################################
###  Fundamental group
################################################################################

@doc raw"""
    fundamental_group(K::SimplicialComplex)

Return the fundamental group of the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> pi_1 = fundamental_group(torus());

julia> describe(pi_1)
"Z x Z"
```
"""
function fundamental_group(K::SimplicialComplex)
    ngens, relations = pm_object(K).FUNDAMENTAL_GROUP
    F = free_group(ngens)
    nrels = length(relations)
    if nrels==0
        return F
    else
        rvec = Vector{FPGroupElem}(undef, nrels)
        for (i, relation) in enumerate(relations)
            relem = one(F)
            for term in relation
                relem *= F[first(term) + 1]^Int(last(term))
            end
            rvec[i] = relem
        end
        pi_1, _ = quo(F, rvec)
        return pi_1
    end
end

################################################################################
###  Sphere recognition heuristics
################################################################################

"""
    is_sphere(K::SimplicialComplex)

Heuristically check if the abstract simplicial complex `K` is a combinatorial sphere; see [JLLT22](@cite).
Note that this is undecidable in general.
Return `true` if recognized as a sphere.
Return `false` if not a sphere.
Return `nothing` if heuristics unsuccessful.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> is_sphere(K)
false
```
"""
is_sphere(K::SimplicialComplex) = pm_object(K).SPHERE::Union{Bool,Nothing}

"""
    is_ball(K::SimplicialComplex)

Heuristically check if the abstract simplicial complex `K` is a combinatorial ball; see [JLLT22](@cite).
Note that this is undecidable in general.
Return `true` if recognized as a ball.
Return `false` if not a ball.
Return `nothing` if heuristics unsuccessful.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> is_ball(K)
true
```
"""
is_ball(K::SimplicialComplex) = pm_object(K).BALL::Union{Bool,Nothing}

"""
    is_manifold(K::SimplicialComplex)

Check if the abstract simplicial complex `K` is a combinatorial manifold, possibly with boundary.
Note that this is undecidable in general.
Return `true` if recognized as a manifold.
Return `false` if not a manifold.
Return `nothing` if heuristics unsuccessful.

# Examples
```jldoctest
julia> is_manifold(torus())
true

julia> is_manifold(simplicial_complex([[1,2],[2,3]]))
true

julia> is_manifold(simplicial_complex([[1,2],[2,3],[2,4]]))
false
```
"""
is_manifold(K::SimplicialComplex) = pm_object(K).MANIFOLD::Union{Bool,Nothing}

################################################################################
###  Surface examples
################################################################################

"""
    torus()

Construct Möbius' (vertex-minimal) 7-vertex triangulation of the torus (surface).
"""
torus() = SimplicialComplex(Polymake.topaz.torus())

"""
    klein_bottle()

Construct a 9-vertex triangulation of the Klein bottle.
"""
klein_bottle() = SimplicialComplex(Polymake.topaz.klein_bottle())

"""
    real_projective_plane()

Construct the (vertex-minimal) 6-vertex triangulation of the real projective plane.
"""
real_projective_plane() = SimplicialComplex(Polymake.topaz.real_projective_plane())

################################################################################
###  Other examples
################################################################################

"""
    complex_projective_plane()

Construct the (vertex-minimal) 9-vertex triangulation of the complex projective plane.
"""
complex_projective_plane() = SimplicialComplex(Polymake.topaz.complex_projective_plane())

################################################################################
###  Subcomplexes
################################################################################

@doc raw"""
    star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the star of the face `sigma` in the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> star_subcomplex(K,[1])
Abstract simplicial complex of dimension 2 on 3 vertices
```
"""
star_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}}) = SimplicialComplex(Polymake.topaz.star_subcomplex(pm_object(K), Polymake.to_zero_based_indexing(sigma)))

@doc raw"""
    link_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the link of the face `sigma` in the abstract simplicial complex `K`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> link_subcomplex(K,[2,3])
Abstract simplicial complex of dimension 0 on 2 vertices
```
"""
link_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}}) = SimplicialComplex(Polymake.topaz.link_subcomplex(pm_object(K), Polymake.to_zero_based_indexing(sigma)))

@doc raw"""
    induced_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})

Return the subcomplex of the abstract simplicial complex `K` consisting of all faces
contained in the nonempty vertex set `sigma`.

The vertices of the result are numbered consecutively from 1, following the ascending order
of `sigma`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1,2,3],[2,3,4]]);

julia> L = induced_subcomplex(K, [1,2,4])
Abstract simplicial complex of dimension 1 on 3 vertices

julia> sort(sort.(collect.(facets(L))))
2-element Vector{Vector{Int64}}:
 [1, 2]
 [2, 3]
```
"""
function induced_subcomplex(K::SimplicialComplex, sigma::Union{Vector{Int}, Set{Int}})
  @req !isempty(sigma) "The vertex set must not be empty"
  @req all(in(1:n_vertices(K)), sigma) "The vertices must be in 1:$(n_vertices(K))"
  vertices = Polymake.to_zero_based_indexing(sigma)
  return SimplicialComplex(Polymake.topaz.induced_subcomplex(pm_object(K), vertices))
end

###############################################################################
### Display
###############################################################################

function Base.show(io::IO, K::SimplicialComplex)
    d = dim(K)
    n = n_vertices(K)
    print(io, "Abstract simplicial complex of dimension $(d) on $(n) vertices")
end

###############################################################################
### Isomorphism
###############################################################################

@doc raw"""
    is_isomorphic(K1::SimplicialComplex, K2::SimplicialComplex)

Check if the given simplicial complexes are isomorphic.

# Examples
```jldoctest
julia> K1 = simplicial_complex([[1,2,3],[2,3,4]]);

julia> K2 = simplicial_complex([[1,2,3],[2,3,4]]);

julia> is_isomorphic(K1, K2)
true
```
"""
function is_isomorphic(K1::SimplicialComplex, K2::SimplicialComplex)
  return Polymake.topaz.isomorphic(pm_object(K1), pm_object(K2))::Bool
end

###############################################################################
### helpful functions
###############################################################################

@doc raw"""
    connected_sum(K1::SimplicialComplex, K2::SimplicialComplex, f1::Int=0, f2::Int=0)

Compute the connected sum of two abstract simplicial complexes. Parameters `f1` and `f2` specify which facet
 of the first and second complex correspondingly are glued together. Default is the
0-th facet of both. The vertices in the selected facets are identified with each
other according to their order in the facet (that is, in increasing index order).

# Examples
```jldoctest
julia> K = torus();

julia> surface_genus_2 = connected_sum(K, K)
Abstract simplicial complex of dimension 2 on 11 vertices

julia> homology(surface_genus_2, 1)
Z^4

julia> is_manifold(surface_genus_2)
true
```
"""
function connected_sum(K1::SimplicialComplex, K2::SimplicialComplex, f1::Int=0, f2::Int=0)
  return SimplicialComplex(Polymake.topaz.connected_sum(pm_object(K1), pm_object(K2), f1, f2))
end

@doc raw"""
    deletion(K::SimplicialComplex, face::Union{<:AbstractSet{Int},<:AbstractVector{Int}})

Remove the given face and all the faces containing it from an abstract simplicial complex `K`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2, 3], [2, 3, 4]]);

julia> K_with_deletion = deletion(K, Set([1, 2]));

julia> facets(K_with_deletion)
2-element Vector{Set{Int64}}:
 Set([3, 1])
 Set([4, 2, 3])
```
"""
function deletion(K::SimplicialComplex, face::Union{<:AbstractSet{Int},<:AbstractVector{Int}})
  zero_based = Polymake.to_zero_based_indexing(face)
  return SimplicialComplex(Polymake.topaz.deletion(pm_object(K), zero_based))
end

@doc raw"""
    automorphism_group(K::SimplicialComplex; action=:on_vertices)

Given a simplicial complex `K` return its automorphism group as a `PermGroup`.
The group can be returned as a subgroup of the permutation group of the vertices
by passing `:on_vertices` to the `action` keyword argument or on the facets
 by passing `:on_facets`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1, 2, 3], [2, 3, 4]])
Abstract simplicial complex of dimension 2 on 4 vertices

julia> automorphism_group(K)
Permutation group of degree 4
```
"""
function automorphism_group(K::SimplicialComplex; action=:on_vertices)
  pm_K = Oscar.pm_object(K)
  Polymake.topaz.combinatorial_symmetries(pm_K)
  if action == :on_vertices
    gens_G = Polymake.to_one_based_indexing(pm_K.GROUP.RAYS_ACTION.GENERATORS)
    n = n_vertices(K)
    return permutation_group(n, perm.(gens_G))
  elseif action == :on_facets
    gens_G = Polymake.to_one_based_indexing(pm_K.GROUP.FACETS_ACTION.GENERATORS)
    n = n_facets(K)
    return permutation_group(n, perm.(gens_G))
  else
    error("unsupported keyword passed to action")
  end
end

@doc raw"""
    on_simplicial_complex(K::SimplicialComplex, g::PermGroupElem)

Given a simplicial complex `K` return the simplicial complex corresponding
to a permutation on it's vertices given by `g`.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2, 3], [2, 3, 4]])
Abstract simplicial complex of dimension 2 on 4 vertices

julia> G = automorphism_group(K)
Permutation group of degree 4

julia> g = collect(G)[2]
(1,4)

julia> facets(on_simplicial_complex(K, g))
2-element Vector{Set{Int64}}:
 Set([2, 3, 1])
 Set([4, 2, 3])
```
"""
function on_simplicial_complex(K::SimplicialComplex, g::PermGroupElem)
  @req degree(parent(g)) == n_vertices(K) "g needs to be an element of the permutation group on the vertices"
  new_facets = on_sets_sets(Set(facets(K)), g)
  simplicial_complex(collect(new_facets))
end

@doc raw"""
     simplicial_product(K1::SimplicialComplex, K2::SimplicialComplex)

Given simplicial complexes `K1` and `K2` return their simplicial product.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K1 = simplicial_complex([[1, 2], [2, 3]])
Abstract simplicial complex of dimension 1 on 3 vertices

julia> K2 = simplicial_complex([[1, 2], [1, 3]])
Abstract simplicial complex of dimension 1 on 3 vertices

julia> facets(simplicial_product(K1, K2))
8-element Vector{Set{Int64}}:
 Set([5, 2, 1])
 Set([5, 4, 1])
 Set([2, 8, 1])
 Set([7, 8, 1])
 Set([6, 2, 3])
 Set([5, 6, 2])
 Set([2, 9, 3])
 Set([2, 9, 8])
```
"""
function simplicial_product(K1::SimplicialComplex, K2::SimplicialComplex)
  return SimplicialComplex(Polymake.topaz.simplicial_product(pm_object(K1), pm_object(K2)))
end

@doc raw"""
     link_subcomplex(K::SimplicialComplex, f::Union{<:AbstractSet{Int},<:AbstractVector{Int}}))

Given simplicial complex `K` and a face `f` of `K` return the link of `f`.

# Examples
```jldoctest
julia> K = simplicial_complex([[1, 2], [2, 3], [3, 4]])
Abstract simplicial complex of dimension 1 on 4 vertices

julia> facets(link_subcomplex(K, [2]))
2-element Vector{Set{Int64}}:
 Set([1])
 Set([2])
```
"""
function link_subcomplex(K::SimplicialComplex, face::Union{<:AbstractSet{Int},<:AbstractVector{Int}})
  zero_based_face = Polymake.to_zero_based_indexing(face)
  return SimplicialComplex(Polymake.link_subcomplex(pm_object(K), zero_based_face))
end

@doc raw"""
   barycentric_subdivision(K::SimplicialComplex)

Given simplicial complex `K` returns its barycentric subdivision.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = simplicial_complex([[1, 2, 3]])
Abstract simplicial complex of dimension 2 on 3 vertices

julia> facets(barycentric_subdivision(K))
6-element Vector{Set{Int64}}:
 Set([5, 2, 1])
 Set([5, 3, 1])
 Set([6, 2, 1])
 Set([4, 6, 1])
 Set([7, 3, 1])
 Set([4, 7, 1])
```
"""
function barycentric_subdivision(K::SimplicialComplex)
  return SimplicialComplex(Polymake.topaz.barycentric_subdivision(pm_object(K))) 
end

function is_shifted(K::SimplicialComplex)
  return pm_object(K).SHIFTED
end
