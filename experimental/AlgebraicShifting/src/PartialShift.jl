################################################################################
#  helpful functions for  Partial Shifts
################################################################################

const ComplexOrHypergraph = Union{UniformHypergraph, SimplicialComplex}

################################################################################
# Matrix construction helper functions
inversions(g::PermGroupElem) = [(i,j) for j in 2:degree(parent(g)) for i in 1:j-1 if g(i) > g(j)]

function generic_unipotent_matrix(R::MPolyRing)
  x = gens(R)
  n_vars = length(x)
  @req is_square(n_vars) "indeterminants should come from a square matrix"
  n = Int(sqrt(n_vars))
  x = matrix(R, reshape(x, n, n))
  u = identity_matrix(R, n)
  # make unipotent matrix
  for i in 1:n
    for j in i + 1:n
      u[i, j] = x[i, j]
    end
  end
  return u
end

@doc raw"""
    generic_unipotent_matrix(R::MPolyRing)
    generic_unipotent_matrix(F::Field, n::Int)

Construct a unipotent matrix with entries in a polynomial ring `R`.
One can also provide a field `F` and an integer `n`,
then the entries of the unipotent matrix will lie in a multivariate
polynomial ring over `F` with `n^2` variables.

# Examples
```jldoctest
julia> R, x = polynomial_ring(QQ, :x=> (1:2, 1:2))
(Multivariate polynomial ring in 4 variables over QQ, QQMPolyRingElem[x[1, 1] x[1, 2]; x[2, 1] x[2, 2]])

julia> generic_unipotent_matrix(R)
[1   x[1, 2]]
[0         1]

julia> generic_unipotent_matrix(GF(2), 2)
[1   x[1, 2]]
[0         1]
```
"""
function generic_unipotent_matrix(F::Field, n::Int)
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n); cached=false)
  return generic_unipotent_matrix(Fx)
end

@doc raw"""
    rothe_matrix(F::Field, w::WeylGroupElem)
    rothe_matrix(F::Field, p::PermGroupElem)
    rothe_matrix(R::MPolyRing, w::WeylGroupElem)
    rothe_matrix(R::MPolyRing, p::PermGroupElem)

For a base field `F` and a Weyl group element `w` return the matrix with entries in the
multivariate polynomial ring `R` with `n^2` many indeterminants where `n - 1` is the rank of the
root system of the Weyl group.
As `general_linear_group(n^2, R)` has a Bruhat decomposition, any element lies in a unique double coset $BwB$, where $B$ is the Borel group of upper triangular matrices.
The Rothe matrix is a normal form for the matrix on the left of a representative for the double coset corresponding to `w`.
There is also the possibility to pass the underlying polynomial ring `R` instead.
This will be explained further once the corresponding preprint is on the arXiv.
We use the name Rothe matrix because of its resemblance with a Rothe diagram, see [Knu98](@cite).

# Examples
```jldoctest
julia> W = weyl_group(:A, 4)
Weyl group
  of root system of rank 4
    of type A4

julia> s = gens(W)
4-element Vector{WeylGroupElem}:
 s1
 s2
 s3
 s4

julia> w = s[2] * s[3] * s[4]
s2 * s3 * s4

julia> rothe_matrix(GF(2), w)
[1         0         0         0   0]
[0   x[2, 3]   x[2, 4]   x[2, 5]   1]
[0         1         0         0   0]
[0         0         1         0   0]
[0         0         0         1   0]

julia> rothe_matrix(QQ, perm([2, 3, 1]))
[x[1, 3]   1   0]
[x[2, 3]   0   1]
[      1   0   0]

julia> Fx, x = polynomial_ring(GF(2), :x => (1:5, 1:5))
(Multivariate polynomial ring in 25 variables over GF(2), FqMPolyRingElem[x[1, 1] x[1, 2] … x[1, 4] x[1, 5]; x[2, 1] x[2, 2] … x[2, 4] x[2, 5]; … ; x[4, 1] x[4, 2] … x[4, 4] x[4, 5]; x[5, 1] x[5, 2] … x[5, 4] x[5, 5]])

julia> rothe_matrix(Fx, w)
[1         0         0         0   0]
[0   x[2, 3]   x[2, 4]   x[2, 5]   1]
[0         1         0         0   0]
[0         0         1         0   0]
[0         0         0         1   0]

julia> rothe_matrix(Fx, perm([1, 3, 2, 5, 4]))
[1         0   0         0   0]
[0   x[2, 3]   1         0   0]
[0         1   0         0   0]
[0         0   0   x[4, 5]   1]
[0         0   0         1   0]
```
"""
function rothe_matrix(F::Field, w::WeylGroupElem; kw...)
  W = parent(w)
  phi = isomorphism(PermGroup, W)
  return rothe_matrix(F, phi(w); kw...)
end

function rothe_matrix(F::Field, p::PermGroupElem; kw...)
  n = degree(parent(p))
  Fx, _ = polynomial_ring(F, :x => (1:n, 1:n); cached=false)
  return rothe_matrix(Fx, p; kw...)
end

function rothe_matrix(R::MPolyRing{T}, w::WeylGroupElem; kw...) where T
  W = parent(w)
  phi = isomorphism(PermGroup, W)
  return rothe_matrix(R, phi(w))
end

function rothe_matrix(R::MPolyRing{T}, p::PermGroupElem;
                      uhg::Union{UniformHypergraph, Nothing}=nothing) where T
  n = degree(parent(p))
  @req ngens(R) == n^2 "The number of generators of the ring should match the square of the degree of the permutation group"
  u = identity_matrix(R, n)
  x = reshape(gens(R), n, n)
  for (i, j) in inversions(p)
    u[i, j] = !isnothing(uhg) && _set_to_zero(uhg, (i, j)) ? zero(R) : x[i, j]
  end
  return u * permutation_matrix(R, p)
end

################################################################################
#TODO we should add bounds so that we dont have to compute determinants we don't need
@doc raw"""
    compound_matrix(m::MatElem, k::Int)
    compound_matrix(p::PermGroupElem, k::Int)
    compound_matrix(w::WeylGroupElem, k::Int)
    compound_matrix(m::MatElem, K::UniformHypergraph)

Given a matrix `m`, return the matrix where each entry is a `k`$\times$`k`-minor of `m`.
The entries of the compound matrix are ordered with respect to the lexicographic order on sets.
When passed a `PermGroupElem` or `WeylGroupElem`, return the compound matrix for their permutation matrix representation.

Alternatively, passing a `UniformHypergraph` `K` will return the compound matrix with entries the `face_size(K)` minors, and restrict the rows to the rows corresponding to `K`.

# Examples
```jldoctest
julia> M = generic_unipotent_matrix(QQ, 3)
[1   x[1, 2]   x[1, 3]]
[0         1   x[2, 3]]
[0         0         1]

julia> compound_matrix(M, 2)
[1   x[2, 3]   x[1, 2]*x[2, 3] - x[1, 3]]
[0         1                     x[1, 2]]
[0         0                           1]

julia> compound_matrix(perm([1, 3, 2]), 2)
[0   1    0]
[1   0    0]
[0   0   -1]

julia> W = weyl_group(:A, 2)
Weyl group
  of root system of rank 2
    of type A2

julia> compound_matrix(longest_element(W), 2)
[ 0    0   -1]
[ 0   -1    0]
[-1    0    0]

julia> K = uniform_hypergraph([[1, 2], [2, 3]])
UniformHypergraph(3, 2, [[1, 2], [2, 3]])

julia> compound_matrix(M, K)
[1   x[2, 3]   x[1, 2]*x[2, 3] - x[1, 3]]
[0         0                           1]
```
"""
function compound_matrix(m::MatElem, K::UniformHypergraph)
  @req size(m,1) == size(m,2) "Only valid for square matrices"
  n = size(m, 1)
  k = face_size(K)
  nCk = sort(subsets(n, k))
  return matrix(base_ring(m), [det(m[row, col]) for row in faces(K), col in nCk])
end

compound_matrix(m::MatElem, k::Int) = compound_matrix(m, uniform_hypergraph(sort(subsets(size(m, 1), k))))
compound_matrix(p::PermGroupElem, k::Int) = compound_matrix(permutation_matrix(ZZ, p), k)

function compound_matrix(w::WeylGroupElem, k::Int)
  iso = isomorphism(PermGroup, parent(w))
  return compound_matrix(permutation_matrix(ZZ, iso(w)), k)
end

# this corresponds to left B invariance
function _set_to_zero(K::SimplicialComplex, indices::Tuple{Int, Int})
  row, col = indices
  row == col && return false
  K_facets = facets(K)
  for facet in K_facets
    # if row index is not in face there is nothing to do
    !(row in facet) && continue
    # replace index from row with index from col in facet
    S = push!(delete!(copy(facet), row), col)
    !any(is_subset(S, check_facet) for check_facet in K_facets) && return false
  end
  return true
end

_set_to_zero(K::UniformHypergraph, indices::Tuple{Int, Int}) = _set_to_zero(simplicial_complex(K), indices)

###############################################################################
# Exterior shift
###############################################################################
function exterior_shift(K::UniformHypergraph, g::MatElem;
                        (ref!)=ModStdQt.ref_ff_rc!,
                        lower_uhg::UniformHypergraph=uniform_hypergraph(Vector{Int}[]))
  # the exterior shifting works in a different algebra that lends
  # itself to an easier implementation
  @req size(g, 1) == size(g, 2) "Change of basis matrix must be square."
  @req size(g, 1) == n_vertices(K) "Matrix size does not match K."
  matrix_base = base_ring(g)
  nCk = collect(combinations(n_vertices(K), face_size(K)))
  k = face_size(K)
  c = compound_matrix(g, K)
  if matrix_base isa MPolyRing
    if !isempty(lower_uhg)
      zero_cols_indices = findall(x -> !all((low_f in faces(lower_uhg)) for low_f in combinations(x, k - 1)), nCk)
      c[:, zero_cols_indices] .= zero(matrix_base)
    end
    ref!(c)
  else
    rref!(c) #TODO could use ref! with different heuristic
  end
  result = uniform_hypergraph(nCk[pivot_columns(c)], n_vertices(K), face_size(K))
  @assert length(faces(result)) == length(faces(K)) "Shifting by invertible `g` should preserve size of hypergraph."
  return result
end

function exterior_shift(K::SimplicialComplex, g::MatElem; kw...)
  n = n_vertices(K)
  f_vec = f_vector(K)
  k = 2
  lower_uhgs = [uniform_hypergraph(Vector{Int}[])]
  while k <= length(f_vec)
    shifted_uhg = exterior_shift(uniform_hypergraph(K, k), g; lower_uhg=lower_uhgs[k - 1], kw...)
    push!(lower_uhgs, shifted_uhg)
    k += 1
  end
  return simplicial_complex([
    (faces(uhg) for uhg in lower_uhgs)...;
    [[i] for i in 1:n] # Make sure result is a complex on n vertices
  ])
end

@doc raw"""
    exterior_shift([F::Field,] K::SimplicialComplex [, w::WeylGroupElem]; las_vegas_trials=100)
    exterior_shift([F::Field,] K::UniformHypergraph [, w::WeylGroupElem]; las_vegas_trials=100)

Compute the (partial) exterior shift of a simplical complex or uniform hypergraph `K` with respect to the Weyl group element `w` and the field `F`.
If the field is not given then `QQ` is used during the computation.
If `w` is not given then `longest_element(weyl_group(:A, n_vertices(K) - 1))` is used.
The algorithm used is a Las Vegas algorithm and you can adjust the number of trials being used with the kwarg `las_vegas_trials` which is set to 100 by default.
Setting `las_vegas_trials=0` will run the deterministic algorithm.

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> K = real_projective_plane()
Abstract simplicial complex of dimension 2 on 6 vertices

julia> is_shifted(K)
false

julia> L = exterior_shift(K)
Abstract simplicial complex of dimension 2 on 6 vertices

julia> facets(L)
10-element Vector{Set{Int64}}:
 Set([2, 3, 1])
 Set([4, 2, 1])
 Set([5, 2, 1])
 Set([6, 2, 1])
 Set([4, 3, 1])
 Set([5, 3, 1])
 Set([6, 3, 1])
 Set([5, 4, 1])
 Set([4, 6, 1])
 Set([5, 6, 1])

julia> is_shifted(L)
true

julia> betti_numbers(L)
3-element Vector{Int64}:
 0
 0
 0

julia> W = weyl_group(:A, n_vertices(K) - 1)
Weyl group
  of root system of rank 5
    of type A5

julia> s = gens(W)
5-element Vector{WeylGroupElem}:
 s1
 s2
 s3
 s4
 s5

julia> w = s[2] * s[3] * s[4]
s2 * s3 * s4

julia> exterior_shift(GF(2), K, w)
Abstract simplicial complex of dimension 2 on 6 vertices
```
"""
function exterior_shift(F::Field, K::ComplexOrHypergraph,
                        p::PermGroupElem; las_vegas_trials::Int=100, kw...)
  n = n_vertices(K)
  @req n == degree(parent(p)) "number of vertices - 1 should equal the rank of the root system"
  las_vegas_trials > 0 && return exterior_shift_lv(F, K, p; n_samples=las_vegas_trials, kw...)
  return exterior_shift(K, rothe_matrix(F, p); kw...)
end

function exterior_shift(F::Field, K::ComplexOrHypergraph, w::WeylGroupElem; kw...)
  phi = isomorphism(PermGroup, parent(w))
  return exterior_shift(F, K, phi(w); kw...)
end

function exterior_shift(F::Field, K::ComplexOrHypergraph; kw...)
  n = n_vertices(K)
  p = perm(reverse(1:n))
  return exterior_shift(F, K, p; kw...)
end

exterior_shift(K::ComplexOrHypergraph; kw...) = exterior_shift(QQ, K; las_vegas_trials=1, kw...)

################################################################################
# Las Vegas Partial Shifting

function random_rothe_matrix(F::Field, p::PermGroupElem)
  n = degree(parent(p))
  u = identity_matrix(F, n)
  for (i, j) in inversions(p)
    u[i, j] = iszero(characteristic(F)) ? F(rand(-10:10)) : rand(F)
  end
  return u * permutation_matrix(F, p)
end

# returns true if the target is the partial shift of src with respect to p
# CAUTION! This function only works correctly if target is obtained as the shift of some matrix.
# this is an improved version of algorithm 3 in dx.doi.org/10.1145/3747199.3747562
function check_shifted(F::Field,
                       src::UniformHypergraph,
                       target::UniformHypergraph,
                       p::PermGroupElem;
                       lower_uhg::UniformHypergraph=uniform_hypergraph(Vector{Int}[]),
                       (ref!)=ModStdQt.ref_ff_rc!)
  # need to check if this sort can be removed
  target_faces = faces(target)
  max_face = length(target_faces) == 1 ? target_faces[1] : max(target_faces...)
  # currently number of faces of src and target are the same
  # this may change in the future
  num_rows = length(faces(src))
  n = n_vertices(src)
  k = face_size(src)
  # check if we can use properties of the target being a shifted family (higher coface elimination)
  full_shift = (p == perm(reverse(1:n)))
  nCk = combinations(n, k)
  # limits the columns by the max face of source
  col_sets = [c for c in nCk if c < max_face]
  # restricted columns is used when we know apriori that certain columns
  # cannot appear in the shifted complex.
  # for example when we know that a column corresponds to a face that contains
  # a lower dimensional non face

  # if the # of non zeros cols up to max_face_index is equal to the rank
  # we do not need to do any row reduction
  n_dependent_columns = length(col_sets) + 1 - num_rows

  if !isempty(lower_uhg)
    zero_cols_indices = findall(x -> !all((low_f in faces(lower_uhg)) for low_f in combinations(x, k - 1)), col_sets)
    needs_check = n_dependent_columns - length(zero_cols_indices) > 0
  else
    zero_cols_indices = Int[]
    needs_check = n_dependent_columns > 0
  end

  if needs_check
    r = rothe_matrix(F, p; uhg=src)
    M = compound_matrix(r, src)[collect(1:num_rows), 1:length(col_sets)]
    if !isempty(zero_cols_indices)
      M[:, zero_cols_indices] .= zero(F)
    end
    
    dep_col_inds = [i for (i, c) in enumerate(col_sets) if !(c in faces(target))]
    cols_to_check = Int[]
    for (i, j) in enumerate(dep_col_inds)
      if !iszero(M[:, j])
        push!(cols_to_check, j)
      end
    end

    if full_shift
      for col in zero_cols_indices
        cols_to_check = [i for i in cols_to_check if col != i && !_domination(col_sets[col], col_sets[i])]
      end
      # set columns we don't need to check to zero
      for i in dep_col_inds
        i in cols_to_check && continue
        M[:, i] .= zero(base_ring(M))
      end
    end
    isempty(cols_to_check) && return true
    max_col = max(cols_to_check...)
    ref!(@view M[:, 1:max_col])
    pivots = pivot_columns(M)
    any([c in pivots for c in cols_to_check]) && return false
  end
  return true
end

# CAUTION! See the comment in the previous functions
function check_shifted(F::Field,
                       src::SimplicialComplex,
                       target::SimplicialComplex,
                       p::PermGroupElem)
  n = n_vertices(src)
  f_vec = f_vector(src)
  k = 2
  lower_uhg = uniform_hypergraph(Vector{Int}[])
  while k <= length(f_vec)
    uhg_src = uniform_hypergraph(src, k)
    uhg_target = uniform_hypergraph(target, k)
    !check_shifted(F, uhg_src, uhg_target, p;
                   lower_uhg=lower_uhg) && return false
    lower_uhg = uhg_target
    k += 1
  end
  return true
end

"""
    exterior_shift_lv(F::Field, K::ComplexOrHypergraph, p::PermGroupElem; n_samples=100)

Computes the (partial) exterior shift of a simplical complex or uniform hypergraph `K` with respect to the permutation group element `p` and the field `F`,
using the Las Vegas algorithm. It samples `n_samples` random matrices, computes the respective partial shifts, takes the lexicographically minimal one,
tests if it is the partial shift of `K` with respect to `p`, and returns the shift.
If the shift was not found using `n_samples`, the algorithm generates `n_samples` in larger and larger field extensions of `F` until found.
"""
function exterior_shift_lv(F::Field, K::ComplexOrHypergraph, p::PermGroupElem; n_samples=100, kw...)
  L = F
  extension = Int(log(characteristic(L), order(L)))
  char = characteristic(F)
  while true
    random_matrices = [random_rothe_matrix(L, p) for _ in 1:n_samples]
    (shift, i) = efindmin((exterior_shift(K, r) for (i, r) in enumerate(random_matrices)); lt=isless_lex) 
    # Check if `shift` is the generic exterior shift of K
    prime_field = iszero(char) ? QQ : fpField(UInt(char))
    n = n_vertices(K)
    is_correct_shift = (p != perm(reverse(1:n)) || is_shifted(shift)) && check_shifted(prime_field, K, shift, p; kw...)

    is_correct_shift && return shift

    # shift not found so extend the field for sampling
    extension += 1
    L = GF(char^extension)
  end
end

function exterior_shift_lv(F::QQField, K::ComplexOrHypergraph, p::PermGroupElem; n_samples=100, kw...)
  while true
    random_matrices = [random_rothe_matrix(F, p) for _ in 1:n_samples]
    (shift, i) = efindmin((exterior_shift(K, r) for (i, r) in enumerate(random_matrices)); lt=isless_lex) 
    # Check if `shift` is the generic exterior shift of K
    n = n_vertices(K)

    is_correct_shift = (p != perm(reverse(1:n)) || is_shifted(shift)) && check_shifted(F, K, shift, p; kw...)

    is_correct_shift && return shift
  end
end
