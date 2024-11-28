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

Constructs a unipotent matrix with entries in a polynomial ring `R`.
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
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
  return generic_unipotent_matrix(Fx)
end

@doc raw"""
    rothe_matrix(F::Field, w::WeylGroupElem; K::Union{SimplicialComplex, Nothing} = nothing)

For a base field `F` and a Weyl group element `w` return the matrix with entries in the
multivariate polynomial ring `R` with `n^2` many indeterminants where `n - 1` is the rank of the
root system of the Weyl group.
As `general_linear_group(n^2, R)` has a Bruhat decomposition, any element lies in a unique double coset $BwB$, where $B$ is the Borel group of upper triangular matrices.
The Rothe matrix is a normal form for the matrix on the left of a representative for the double coset corresponding to `w`.
This will be explained further once the corresponding preprint is on the arXiv.
We use the name Rothe matrix because of its resemblance with a Rothe diagram. (add ref? Knuth?)

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
```
"""
function rothe_matrix(F::Field, w::WeylGroupElem)
  W = parent(w)
  phi = isomorphism(PermGroup, W)
  return rothe_matrix(F, phi(w))
end

function rothe_matrix(F::Field, p::PermGroupElem)
  n = degree(parent(p))
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
  u = identity_matrix(Fx, n)
  for (i, j) in inversions(p)
    u[i, j] = x[i, j]
  end
  return u * permutation_matrix(F, p)
end

@doc raw"""
    compound_matrix(m::MatElem, k::Int)
    compound_matrix(p::PermGroupElem, k::Int)
    compound_matrix(w::WeylGroupElem, k::Int)
    compound_matrix(m::MatElem, K::Vector{Vector{Int}})

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
    exterior_shift(F::Field, K::SimplicialComplex, w::WeylGroupElem)
    exterior_shift(F::Field, K::UniformHypergraph, w::WeylGroupElem)
    exterior_shift(K::SimplicialComplex, w::WeylGroupElem)
    exterior_shift(K::UniformHypergraph, w::WeylGroupElem)
    exterior_shift(K::SimplicialComplex)
    exterior_shift(K::UniformHypergraph)

Computes the (partial) exterior shift of a simplical complex or uniform hypergraph `K` with respect to the Weyl group element `w` and the field `F`.
If the field is not given then `QQ` is used during the computation.
If `w` is not given then `longest_element(weyl_group(:A, n_vertices(K) - 1))` is used

# Examples
```jldoctest
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
function exterior_shift(F::Field, K::ComplexOrHypergraph, p::PermGroupElem)
  n = n_vertices(K)
  @req n == degree(parent(p)) "number of vertices - 1 should equal the rank of the root system"
  # might want to reintroduce these lines at some point
  #bool_mat = matrix(F, [!Oscar._set_to_zero(K, (i, j)) for i in 1:n, j in 1:n])
  #M = (bool_mat * permutation_matrix(F, p)) .* rothe_matrix(F, p)
  return exterior_shift(K, rothe_matrix(F, p))
end

function exterior_shift(F::Field, K::ComplexOrHypergraph, w::WeylGroupElem)
  n = n_vertices(K)
  phi = isomorphism(PermGroup, parent(w))
  return exterior_shift(F, K, phi(w))
end

function exterior_shift(F::Field, K::ComplexOrHypergraph)
  n = n_vertices(K)
  W = weyl_group(:A, n - 1)
  return exterior_shift(F, K, longest_element(W))
end

exterior_shift(K::ComplexOrHypergraph) = exterior_shift(QQ, K)

################################################################################
# Las Vegas Partial Shifting

# returns a random invertible matrices of the given sample size
function random_invertible_matrices(n::Int; range::Int = 100, sample_size::Int=100, F::Field=QQ)
  n_samples = 0
  samples = Set{MatElem}()
  n_non_inv = 0
  while (n_samples < sample_size)
    a = matrix(F, round.(Integer, range .* rand(n,n)) .- range // 2)
    if iszero(det(a))
      n_non_inv += 1
      continue
    end
    push!(samples, a)
    n_samples = length(samples)
  end
  # @info "probability of non invertible matrix" n_non_inv / (n_non_inv + sample_size)
  return collect(samples)
end

function random_unipotent_matrix(F::Field, n::Int)
  char = characteristic(F)
  range = char == 0 ? 100 : char
  upper_triangular_matrix(
    reduce(vcat, [
      [F(1); F.(round.(Integer, range .* rand(i)) .- (range // 2))] for i in reverse(0:n-1)
        ]))
end

function random_shift(F::Field, K::ComplexOrHypergraph, p::PermGroupElem;)
  n = n_vertices(K)
  exterior_shift(K, random_unipotent_matrix(F, n) * permutation_matrix(F, p))
end

random_shift(K::ComplexOrHypergraph, p::PermGroupElem) = random_shift(QQ, K, p)

# returns true if the dst is the partial shift of src with respect to w
function check_shifted(F::Field, w::WeylGroupElem,
                       src::UniformHypergraph, dst::UniformHypergraph)
  dst_faces = faces(dst)
  if length(dst_faces) == 1
    max_face = dst_faces[1]
  else
    max_face = max(dst_faces...)
  end
  num_rows = length(dst_faces)
  n = n_vertices(src)
  k = face_size(src)
  nCk = sort(subsets(n, k))
  max_face_index = findfirst(x -> x == max_face, nCk)
  cols = nCk[1:max_face_index - 1]
  r = rothe_matrix(F, w)
  
  if max_face_index > num_rows 
    M = compound_matrix(r, src)[collect(1:num_rows), collect(1:length(cols))]
    Oscar.ModStdQt.ref_ff_rc!(M)
    nCk[independent_columns(M)] == dst_faces && return false
  end
  return true
end

function check_shifted(F::Field, w::WeylGroupElem,
                       src::SimplicialComplex, dst::SimplicialComplex)
  n = n_vertices(src)
  f_vec = f_vector(src)
  k = length(f_vec)

  while k > 1
    uhg_src = uniform_hypergraph(complex_faces(src, k - 1), n)
    uhg_dst = uniform_hypergraph(complex_faces(dst, k - 1), n)
    !check_shifted(F, w, uhg_src, uhg_dst) && return false
    k -= 1
  end
  return true
end

function partial_ext_shift_lv(F::Field, w::WeylGroupElem, K::ComplexOrHypergraph)
  sample_size = characteristic(F) == 0 ? 10 : 100
  W = parent(w)
  phi = isomorphism(PermGroup, W)

  shift = partialsort!([random_shift(F, K, phi(w)) for _ in 1:sample_size], 1;
                      lt=isless_lex)
  check_shifted(F, w, K, shift) && return shift
  return nothing # partial_ext_shift_lv(F, w, K)
end

