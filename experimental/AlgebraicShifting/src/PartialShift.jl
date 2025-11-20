################################################################################
#  helpful functions for  Partial Shifts
################################################################################

const ComplexOrHypergraph = Union{UniformHypergraph, SimplicialComplex}

################################################################################
# Matrix construction helper functions
function inversions(g::PermGroupElem)
  return [(i,j) for j in 2:degree(parent(g)) for i in 1:j-1 if g(i) > g(j)]
end

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
function rothe_matrix(F::Field, w::WeylGroupElem)
  W = parent(w)
  phi = isomorphism(PermGroup, W)
  return rothe_matrix(F, phi(w))
end

function rothe_matrix(F::Field, p::PermGroupElem)
  n = degree(parent(p))
  Fx, _ = polynomial_ring(F, :x => (1:n, 1:n); cached=false)
  return rothe_matrix(Fx, p)
end

function rothe_matrix(R::MPolyRing{T}, w::WeylGroupElem) where T
  W = parent(w)
  phi = isomorphism(PermGroup, W)
  return rothe_matrix(R, phi(w))
end

function rothe_matrix(R::MPolyRing{T}, p::PermGroupElem) where T
  n = degree(parent(p))
  @req ngens(R) == n^2 "The number of generators of the ring should match the square of the degree of the permutation group"
  u = identity_matrix(R, n)
  x = reshape(gens(R), n, n)
  for (i, j) in inversions(p)
    u[i, j] = x[i, j]
  end
  return u * permutation_matrix(R, p)
end

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

# this might be removed (currently is not used)
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
function exterior_shift(K::UniformHypergraph, g::MatElem)
  # the exterior shifting works in a different algebra that lends
  # itself to an easier implementation
  @req size(g, 1) == size(g, 2) "Change of basis matrix must be square."
  @req size(g, 1) == n_vertices(K) "Matrix size does not match K."
  matrix_base = base_ring(g)
  nCk = sort!(subsets(n_vertices(K), face_size(K)))
  c = compound_matrix(g, K)
  if matrix_base isa MPolyRing
    Oscar.ModStdQt.ref_ff_rc!(c)
  elseif matrix_base isa MPolyQuoRing
    lifted_c = lift.(c)
    Oscar.ModStdQt.ref_ff_rc!(lifted_c)
    c = simplify.(matrix_base.(lifted_c))
  else
    rref!(c)
  end
  return uniform_hypergraph(nCk[independent_columns(c)], n_vertices(K), face_size(K))
end

function exterior_shift(K::SimplicialComplex, g::MatElem)
  return simplicial_complex([
    (faces(exterior_shift(uniform_hypergraph(K, k), g)) for k in 1:dim(K)+1)...;
    [[i] for i in 1:n_vertices(K)] # Make sure result is a complex on n vertices
  ])
end

@doc raw"""
    exterior_shift(F::Field, K::SimplicialComplex, w::WeylGroupElem; las_vegas=false)
    exterior_shift(F::Field, K::UniformHypergraph, w::WeylGroupElem; las_vegas=false)
    exterior_shift(K::SimplicialComplex, w::WeylGroupElem; las_vegas=false)
    exterior_shift(K::UniformHypergraph, w::WeylGroupElem; las_vegas=false)
    exterior_shift(K::SimplicialComplex; las_vegas=false)
    exterior_shift(K::UniformHypergraph; las_vegas=false)

Compute the (partial) exterior shift of a simplical complex or uniform hypergraph `K` with respect to the Weyl group element `w` and the field `F`.
If the field is not given then `QQ` is used during the computation.
If `w` is not given then `longest_element(weyl_group(:A, n_vertices(K) - 1))` is used.
Setting `las_vegas=true` will run the algorithm with a random change of basis matrix and repeat the algorithm until the shift is found.

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

julia> betti_numbers(L) == betti_numbers(K)
true

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

julia> L = exterior_shift(GF(2), K, w)
Abstract simplicial complex of dimension 2 on 6 vertices
```
"""
function exterior_shift(F::Field, K::ComplexOrHypergraph,
                        p::PermGroupElem; las_vegas::Bool=false)
  n = n_vertices(K)
  @req n == degree(parent(p)) "number of vertices - 1 should equal the rank of the root system"
  las_vegas && return exterior_shift_lv(F, K, p)
  return exterior_shift(K, rothe_matrix(F, p))
end

function exterior_shift(F::Field, K::ComplexOrHypergraph, w::WeylGroupElem; kw...)
  n = n_vertices(K)
  phi = isomorphism(PermGroup, parent(w))
  return exterior_shift(F, K, phi(w); kw...)
end

function exterior_shift(F::Field, K::ComplexOrHypergraph; kw...)
  n = n_vertices(K)
  W = weyl_group(:A, n - 1)
  return exterior_shift(F, K, longest_element(W); kw...)
end

exterior_shift(K::ComplexOrHypergraph; kw...) = exterior_shift(QQ, K; kw...)

################################################################################
# Las Vegas Partial Shifting

function random_rothe_matrix(F::QQField, p::PermGroupElem)
  n = degree(parent(p))
  u = identity_matrix(F, n)
  for (i, j) in inversions(p)
    u[i, j] = F(Rational(rand() - 0.5))
  end
  return u * permutation_matrix(F, p)
end

function random_rothe_matrix(F::Field, p::PermGroupElem)
  @req !iszero(characteristic(F)) "Field should have positive characteristic"
  range = characteristic(F)
  n = degree(parent(p))
  u = identity_matrix(F, n)
  for (i, j) in inversions(p)
    u[i, j] = F(rand(1:range))
  end
  return u * permutation_matrix(F, p)
end

function random_shift(F::QQField, K::ComplexOrHypergraph, p::PermGroupElem)
  n = n_vertices(K)
  exterior_shift(K, random_rothe_matrix(F, p))
end

function random_shift(F::Field, K::ComplexOrHypergraph, p::PermGroupElem)
  n = n_vertices(K)
  exterior_shift(K, random_rothe_matrix(F, p))
end

random_shift(K::ComplexOrHypergraph, p::PermGroupElem) = random_shift(QQ, K, p)

# returns true if the target is the partial shift of src with respect to p
function check_shifted(F::Field, src::UniformHypergraph,
                       target::UniformHypergraph, p::PermGroupElem)
  target_faces = sort(faces(target))
  max_face = length(target_faces) == 1 ? target_faces[1] : max(target_faces...)
  # currently number of faces of src and target are the same
  # this may change in the future
  num_rows = length(faces(src)) 
  n = n_vertices(src)
  k = face_size(src)
  nCk = sort(subsets(n, k))
  max_face_index = findfirst(x -> x == max_face, nCk)
  # limits the columns by the max face of source
  cols = nCk[1:max_face_index - 1]
  r = rothe_matrix(F, p)

  if max_face_index > num_rows
    M = compound_matrix(r, src)[collect(1:num_rows), collect(1:length(cols))]
    Oscar.ModStdQt.ref_ff_rc!(M)
    nCk[independent_columns(M)] != target_faces[1:end - 1] && return false
  end
  return true
end

function check_shifted(F::Field, src::SimplicialComplex,
                       target::SimplicialComplex, p::PermGroupElem)
  n = n_vertices(src)
  f_vec = f_vector(src)
  k = length(f_vec)
  while k > 1
    uhg_src = uniform_hypergraph(complex_faces(src, k - 1), n)
    uhg_target = uniform_hypergraph(complex_faces(target, k - 1), n)
    !check_shifted(F, uhg_src, uhg_target, p) && return false
    k -= 1
  end
  return true
end

function exterior_shift_lv(F::Field, K::ComplexOrHypergraph, p::PermGroupElem)
  # this might need to be changed based on the characteristic
  # we expect that the larger the characteristic the smaller the sample needs to be
  # setting to 100 now for good measure
  sample_size = 100
  shift = partialsort!([random_shift(F, K, p) for _ in 1:sample_size], 1;
                       lt=isless_lex)

  check_shifted(F, K, shift, p) && return shift

  # this should be updated to not throw an error
  error("Could not find the full shift using $sample_size samples")
end

function exterior_shift_lv(F::QQField, K::ComplexOrHypergraph, p::PermGroupElem)
  shift = random_shift(F, K, p) 
  count = 1
  while !check_shifted(F, K, shift, p)
    count += 1
    shift = random_shift(F, K, p) 
  end
  @vprint :AlgebraicShifting 1 "Number of random shifts computed: $count\n"
  return shift
end
