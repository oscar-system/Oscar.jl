import AbstractAlgebra.Field
import AbstractAlgebra.Generic.MatSpaceElem


export AlgGroup, subgroup

@Markdown.doc """
    AlgGroup{FieldType, RingType, RingElemType, MatrixType}

An algebraic group G defined as an algebraic subvariety of the space GL(t,𝕜). 
The information stored consists of 
  * the base field 𝕜 of type `FieldType`;
  * a polynomial ring `R` = 𝕜[x₁,…,xₙ] of type `RingType`, 
    the ambient polynomial ring for G seen as an affine variety;
  * an ideal `I` ⊂ `R` of type `MPolyIdeal{RingElemType}` defining G 
    as a subscheme of 𝔸ⁿ(𝕜);
  * a t×t-matrix `Z` of type `MatrixType` with entries in `R` embedding 
    G as an algebraic subscheme of GL(t,𝕜) and endowing it with its group structure.
"""
mutable struct AlgGroup{FieldType, RingType, RingElemType, MatrixType}
  k::FieldType
  R::RingType 			# The underlying polynomial ring
  I::MPolyIdeal{RingElemType}	# The ideal describing the closed subgroup
  Z::MatrixType			# The matrix coordinates for the representation

  function AlgGroup( k::FieldType, R::RingType, I::MPolyIdeal{RingElemType}, Z) where{FieldType<:Field, RingType<:Ring, RingElemType}
    @assert typeof(Z) <: dense_matrix_type(R) "Matrix is not defined over the given ring"
    @assert base_ring(R) === k
    @assert typeof(Z) == dense_matrix_type(R)
    return new{typeof(k), typeof(R), elem_type(R), dense_matrix_type(R)}(k, R, I, Z)
  end
end

AlgGroup(R::RingType, I::MPolyIdeal{RingElemType}, Z) where{RingType<:Ring, RingElemType} = AlgGroup(base_ring(R), R, I, Z)

function AlgGroup(I::MPolyIdeal{RingElemType}, Z) where{RingElemType}
  @assert typeof(Z) == dense_matrix_type(base_ring(I)) 
  @assert base_ring(Z) == base_ring(Z) "matrix is not defined over the same ring as the given ideal"
  return AlgGroup(base_ring(base_ring(Z)), base_ring(Z), I, Z)
end

@Markdown.doc """
    base_field(G::AlgGroup)

Returns the field 𝕜 over which G is defined.
"""
function base_field(G::AlgGroup{FieldType, RingType, RingElemType, MatrixType}) where{
	FieldType<:Field, RingType<:Ring, RingElemType, MatrixType}
  return G.k::FieldType
end

@Markdown.doc """
    ambient_ring(G::AlgGroup)

Returns the polynomial ring `R` used to define `G` as an affine variety.
"""
function ambient_ring(G::AlgGroup{FieldType, RingType, RingElemType, MatrixType}) where{
	FieldType<:Field, RingType<:Ring, RingElemType, MatrixType}
  return G.R::RingType
end

@Markdown.doc """
    defining_ideal(G::AlgGroup)

Returns the ideal `I` used to define `G` as an affine subvariety of 𝔸ⁿ(𝕜).
"""
function defining_ideal(G::AlgGroup{FieldType, RingType, RingElemType, MatrixType}) where{
	FieldType<:Field, RingType<:Ring, RingElemType, MatrixType}
  return G.I::MPolyIdeal{elem_type(G.R)}
end

@Markdown.doc """
    presentation_matrix(G::AlgGroup)

Returns the t×t-matrix `Z` with entries in `R` for the variety G = Spec `R`/`I` 
which presents G as a subgroup of GL(t,𝕜).
"""
function presentation_matrix(G::AlgGroup{FieldType, RingType, RingElemType, MatrixType}) where{
	FieldType<:Field, RingType<:Ring, RingElemType, MatrixType}
  return G.Z::dense_matrix_type(G.R)
end

@Markdown.doc """
    function GL( n::Int, k<:Field )

The general linear group of invertible n×n-matrices over the field `k`.
"""
function GL( n::Int, k::FieldType ) where{FieldType<:Field}
  #R_t,t = PolynomialRing(QQ,"t")
  R, z,t = PolynomialRing(k, "a" => (1:n, 1:n), "t" => (1:1))
  Z = matrix(z)
  t = t[1]
  I = ideal(R, [1-det(Z)*t] )
  return AlgGroup{typeof(k), typeof(R), elem_type(R), dense_matrix_type(R)}( k, R, I, Z )
end


function subgroup( G::AlgGroup, I::MPolyIdeal )
  @assert parent(I) == parent(G.R) "given ideal does not live in the same ring as the defining ideal of the algebraic group"
  return AlgGroup( G.R, G.I + I, G.Z )
end

@Markdown.doc """
    LinearRepresentation

A linear representation of an algebraic group G over a field 𝕜, 
defined on a vector space V ≅ 𝕜ʳ. 
The information stored consists of 
  * an instance `G` of AlgGroup which describes the group as a subscheme of GL(t,𝕜);
  * an r×r-matrix A with entries in the polynomial ring 𝕜[zᵢⱼ| 1 ≤ i,j ≤ t] 
    inducing the map G →  GL(r,𝕜) ≅ End(V) for the representation.
"""
mutable struct LinearRepresentation{FieldType, RingType, RingElemType, MatrixTypeG, MatrixTypeA}
  G::AlgGroup{FieldType, RingType, RingElemType, MatrixTypeG}
  A::MatrixTypeA

  function LinearRepresentation( G::AlgGroup, A::MatSpaceElem )
    @assert G.k === base_ring(base_ring(A))
    @assert ncols(G.Z)*nrows(G.Z) == length(gens(base_ring(A)))
    return new{typeof(G.k), typeof(G.R), elem_type(G.R), typeof(G.Z), typeof(A)}( G, A )
  end

end

# the identity representation
function LinearRepresentation( G::AlgGroup )
  n = ncols(G.Z)
  S, a = PolynomialRing(base_field(G), "a" => (1:n, 1:n))
  A = matrix(a)
  return LinearRepresentation( G, A )
end

function representation_on_sym_power( rep::LinearRepresentation, d::Int )
  # return the representation induced on the d-th symmetric power 
  # of the vector space V on which the representation rep is 
  # defined.
  n = ncols( rep.A )	# the rank of the original representation
  N = binomial( n+d-1, d )	# the rank of the induced representation
  # A polynomial ring containing variables for the coordinates of 
  # the original representation:
  S, y = PolynomialRing( rep.G.R, [ "y$k" for k in (1:n) ])
  # An empty matrix for the induced representation
  M = MatrixSpace( rep.G.R, N, N )
  B = zero(M)
  # The matrix of the original representation, but promoted to the 
  # new ring S
  A = matrix( S.(rep.A) )
  # images of the variables under the action of the group
  z = A*matrix(y)
  # Populate the matrix B for the induced representation...
  for i in HomogMultiindex( n, d )
    # ...by going through the multiindices of the monomials 
    # in the given degree and taking the respective powers...
    f = power( [z[k,1] for k in (1:n)], i )
    for j in HomogMultiindex( n, d )
      # ...and extracting the coefficients again.
      B[ linear_index(j), linear_index(i) ] = coeff( f, power( y, j ))
    end
  end
  return LinearRep( rep.G, B )
end
