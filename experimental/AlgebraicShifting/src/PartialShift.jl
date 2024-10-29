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

function generic_unipotent_matrix(F::Field, n::Int)
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
  return generic_unipotent_matrix(Fx)
end

@doc raw"""
     generic_unipotent_matrix(R::MPolyRing)
     generic_unipotent_matrix(F::Field, n::Int)

Constructs a unipotent matrix with entries in a polynomial ring `R`.
One can also provide a field `F` and an integer `n`,
then the entries of the unipotent matrix will lie in a multi variate
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

@doc raw"""
     rothe_matrix(F::Field, w::WeylGroupElem; K::Union{SimplicialComplex, Nothing} = nothing)

For a base field `F` and a weyl group element `w` return the matrix with entries in the
multivariate polynomial ring `R` with `n^2` many indeterminants where `n - 1` is the rank of the
root system of the weyl group.
We know that since `general_linear_group(n^2, R)` has a Bruhat decomposition, and element lies in some double coset $BwB$.
The \emph{Rothe matrix} is a normal form for a the matrix on the left of a representative for the double coset corresponding to `w`.
We use the name \emph{Rothe matrix} because of its resemblance with a \emph{Rothe diagram} (add ref?)

# Examples
```jldoctest

```
"""

function rothe_matrix(F::Field, w::WeylGroupElem)
  n = rank(root_system(parent(w)))+1
  Fx, x = polynomial_ring(F, :x => (1:n, 1:n))
  u = identity_matrix(Fx, n)
  for (i, j) in inversions(perm(w))
    u[i, j] = x[i, j]
  end
  return u * permutation_matrix(F, perm(w))
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

""" Given `K` checks if matrix entry can be set to zero zero """
function _set_to_zero(K::SimplicialComplex, indices::Tuple{Int, Int})
  row, col = indices
  row == col && return false
  K_facets = facets(K)
  for facet in K_facets
    # if column index is not in face there is nothing to do
    !(col in facet) && continue
    # replace index from row with index from col in facet
    S = push!(delete!(copy(facet), col), row)
    !any(is_subset(S, check_facet) for check_facet in K_facets) && return false
  end
  return true
end

_set_to_zero(K::UniformHypergraph, indices::Tuple{Int, Int}) = _set_to_zero(simplicial_complex(K), indices)

###############################################################################
# Exterior shift 
###############################################################################
raw"""    exterior_shift(K::UniformHypergraph, g::MatElem)
    
Computes the exterior shift ``\Delta_g(K)`` of ``K`` w.r.t. the invertible matrix ``g``.
"""
function exterior_shift(K::UniformHypergraph, g::MatElem)
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
    (faces(exterior_shift(uniform_hypergraph(K, k), g)) for k in 1:dim(K)+1)...;
    [[i] for i in 1:n_vertices(K)] # Make sure result is a complex on n vertices
  ])
end
  
@doc raw"""
     exterior_shift(F::Field, K::ComplexOrHypergraph, w::WeylGroupElem)
     exterior_shift(K::ComplexOrHypergraph, w::WeylGroupElem)
     exterior_shift(K::ComplexOrHypergraph)

Computes the (partial) exterior shift of a simplical complex or uniform hypergraph `K` with respect to the Weyl group element `w` and the field `F`.
If the field is not given then `QQ` is used during the computation.
If `w` is not given then `longest_element(weyl_group(:A, n_vertices(K) - 1))` is used


# Example
Compute the exterior generic shift of the real projective plane:
```
julia> K = real_projective_plane()
julia> is_shifted(K)
julia> L = exterior_shift(K)
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

  # in some cases this provides a speed up, not extremely important, we can chose to omit
  # these lines
  # handle some pre computation row reduction of compound matrix
  #bool_mat = matrix(F, [!_set_to_zero(K, (i, j)) for i in 1:n, j in 1:n])
  #M = (bool_mat * permutation_matrix(F, perm(w))) .* 
  
  return exterior_shift(K, rothe_matrix(F, w))
end

function exterior_shift(F::Field, K::ComplexOrHypergraph)
  n = n_vertices(K)
  W = weyl_group(:A, n - 1)
  return exterior_shift(F, K, longest_element(W))
end

exterior_shift(K::ComplexOrHypergraph) = exterior_shift(QQ, K)
