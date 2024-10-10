################################################################################
#  helpful functions for  Partial Shifts
################################################################################

const ComplexOrHypergraph = Union{UniformHypergraph, SimplicialComplex}

################################################################################
# Matrix construction helper functions
function inversions(g::PermGroupElem)
  return [(i,j) for j in 2:degree(parent(g)) for i in 1:j-1 if g(i) > g(j)]
end

""" Given `K` checks if index can be set to zero (B- invariance)"""
function is_zero_entry(K::SimplicialComplex, indices::Tuple{Int, Int})
  K_facets = Set{Set{Int}}(facets(K))
  for facet in K_facets
    # need a comment for this line
    !(indices[1] in facet) && continue
    # swap indices[1] with indices[2] in facet
    S = push!(delete!(copy(facet), indices[1]), indices[2])
    !any(is_subset(S, check_facet) for check_facet in K_facets) && return false
  end
  return true
end

is_zero_entry(K::UniformHypergraph, indices::Tuple{Int, Int}) = is_zero_entry(simplicial_complex(K), indices)

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

function generic_unipotent_matrix(F::Field, n::Int)
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
  return generic_unipotent_matrix(Fx)
end

raw""" Returns the matrix ``\frU(w)``. """
function generic_unipotent_matrix(F::Field, w::WeylGroupElem)
  n = rank(root_system(parent(w)))+1
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
  u = identity_matrix(Fx, n)
  for (i, j) in inversions(perm(w))
    u[i, j] = x[i, j]
  end
  return u
end

""" The `K`-indexed rows of the `k`th compound matrix of a square matrix `X`, where `K` is `k`-homogeneous. """
function compound_matrix(m::MatElem, K::Vector{Vector{Int}})	
	@req size(m,1) == size(m,2) "Only valid for square matrices"
	@req length(Set(map(length, K))) == 1 "All entries in K must have the same size."
	n = size(m, 1)
	k = collect(Set(map(length, K)))[1]
	@req all(1 <= i <= n for s in K for i in s) "All entries in K must represent $k-element subsets of [n]."
  nCk = sort(subsets(n, k))
	return matrix(base_ring(m), [det(m[row, col]) for row in K, col in nCk])
end

compound_matrix(m::MatElem, k::Int) = compound_matrix(m, sort(subsets(size(m, 1), k)))
compound_matrix(p::PermGroupElem, k::Int) = compound_matrix(permutation_matrix(ZZ, p), k)
compound_matrix(w::WeylGroupElem, k::Int) = compound_matrix(perm(w), k)

###############################################################################
# Exterior shift 
###############################################################################
raw"""    exterior_shift(K::UniformHypergraph, g::MatElem)
    
Computes the exterior shift ``\Delta_g(K)`` of ``K`` w.r.t. the invertible matrix ``g``.
"""
function exterior_shift(K::UniformHypergraph, g::MatElem{})
  # the exterior shifting works in a different algebra that lends
  # itself to an easier implementation 
  @req size(g, 1) == size(g, 2) "Change of basis matrix must be square."
  @req size(g, 1) == n_vertices(K) "Matrix size does not match K."
  matrix_base = base_ring(g)
  nCk = sort!(subsets(n_vertices(K), face_size(K)))
  c = compound_matrix(g, faces(K))
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
    (exterior_shift(faces(uniform_hypergraph(K, k), g)) for k in 1:dim(K)+1)...;
    [[i] for i in 1:n_vertices(K)] # Make sure result is a complex on n vertices
  ])
end
  
raw"""    exterior_shift(K::ComplexOrHypergraph, w :: WeylGroupElem  = nothing;  field :: Field = QQ)

Computes the (partial) generic exterior shift of ``K`` w.r.t. the Weyl group element ``w``. 

With ``w = w_0`` the longest word in the Weyl group (default), this is the generic exterior shift of ``K``.

# Example
Compute the exterior generic shift of the real projective plane:
```
julia> K = real_projective_plane()
julia> is_shifted(K)
julia> L = exterior_shift2(K)
julia> is_shifted(L)
julia> betti_numbers(K) == betti_numbers(L)
```
Apply the partial generic shift w.r.t. a permutation ``w``:
```
julia> W = weyl_group(:A, 5)
julia> s = gens(W)
julia> w = s[1] * s[2] * s[1]
julia> L = exterior_shift(K, w)
```
"""
function exterior_shift(F::Field, K::ComplexOrHypergraph, w::WeylGroupElem)
  n = n_vertices(K)
  @req n == rank(root_system(parent(w))) + 1 "number of vertices - 1 should equal the rank of the root system"

  # set certain entries to zero, see B-invariance (probably needs a better name)
  bool_mat = matrix(F, [is_zero_entry(K, (i, j) for i in 1:n, j in 1:n)])
  M = bool_mat .* generic_unipotent_matrix(F, w)
  return exterior_shift(K,  M * permutation_matrix(F, perm(w)))
end

function exterior_shift(F::Field, K::ComplexOrHypergraph)
  n = n_vertices(K)
  W = weyl_group(:A, n - 1)
  return exterior_shift(F, K, longest_element(W))
end
