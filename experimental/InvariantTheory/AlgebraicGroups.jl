import AbstractAlgebra.Field
import AbstractAlgebra.Generic.MatSpaceElem


export AlgGroup, subgroup, base_field, ambient_ring, defining_ideal, presentation_matrix, AlgGroup_GL, AlgGroup_SL

export LinearRepresentation, representation_on_sym_power, mapping_matrix, group


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

  function AlgGroup( k::FieldType, R::RingType, I::MPolyIdeal{RingElemType}, Z::MatrixType) where{FieldType<:Field, RingType<:Ring, RingElemType, MatrixType}
    @assert typeof(Z) <: dense_matrix_type(R) "Matrix is not defined over the given ring"
    @assert base_ring(R) === k "ambient polynomial ring is not defined over the given field" 
    @assert MatrixType == dense_matrix_type(R) "Matrix has the wrong type"
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
    GL( n::Int, k::FieldType ) where{FieldType<:Field}

The general linear group of invertible n×n-matrices over the field `k`.
"""
function AlgGroup_GL( n::Int, k::FieldType ) where{FieldType<:Field}
  #R_t,t = PolynomialRing(QQ,"t")
  R, z,t = PolynomialRing(k, "z" => (1:n, 1:n), "t" => (1:1))
  Z = matrix(z)
  t = t[1]
  I = ideal(R, [1-det(Z)*t] )
  return AlgGroup( k, R, I, Z )
end

@Markdown.doc """
    function SL(n::Int, k::FieldType) where{FieldType<:Field}

Constructs the special linear group of rank `n` over the field `k`.
"""
function AlgGroup_SL(n::Int, k::FieldType) where{FieldType<:Field}
  R, z = PolynomialRing(k, "z" => (1:n, 1:n))
  Z = matrix(z)
  I = ideal(R, [1-det(Z)] )
  return AlgGroup( k, R, I, Z )
end
  


function subgroup( G::AlgGroup, I::MPolyIdeal )
  @assert parent(I) == parent(G.R) "given ideal does not live in the same ring as the defining ideal of the algebraic group"
  return AlgGroup( G.R, G.I + I, G.Z )
end

@Markdown.doc """
    LinearRepresentation{FieldType, RingType, RingElemType, MatrixTypeG, MatrixTypeA}

A linear representation of an algebraic group G ⊂ GL(t,𝕜) over a field 𝕜, 
defined on a vector space V ≅ 𝕜ʳ. 
The information stored consists of 
  * an instance `G` of `AlgGroup{FieldType, RingType, RingElemType, MatrixTypeG}` 
    which describes the group as a subscheme of GL(t,𝕜);
  * an r×r-matrix A of type `MatrixTypeA` with entries in the polynomial ring 𝕜[zᵢⱼ| 1 ≤ i,j ≤ t],
    whose variables are the standard coordinates on GL(t,𝕜), inducing the 
    map G →  GL(r,𝕜) ≅ End(V) for the representation.
"""
mutable struct LinearRepresentation{FieldType, RingType, RingElemType, MatrixTypeG, MatrixTypeA}
  G::AlgGroup{FieldType, RingType, RingElemType, MatrixTypeG}
  A::MatrixTypeA
  pullback_morphism::Oscar.AlgHom
  #pullback_morphism::AlgebraHomomorphism

  function LinearRepresentation( G::AlgGroup, A::MatSpaceElem )
    @assert G.k === base_ring(base_ring(A))
    @assert ncols(G.Z)*nrows(G.Z) == length(gens(base_ring(A)))
    images = Vector{elem_type(ambient_ring(G))}()
    Z = presentation_matrix(G)
    n = ncols(Z)
    # calling gens(base_ring(A)) reveals that the order of the 
    # coordinates is reversed...
    for j in (1:n)
      for i in (1:n)
	push!(images, Z[i,j])
      end
    end
    phi = AlgebraHomomorphism(base_ring(A), ambient_ring(G), images)
    return new{typeof(G.k), typeof(G.R), elem_type(G.R), typeof(G.Z), typeof(A)}(G, A, phi)
  end

end

# the identity representation
function LinearRepresentation( G::AlgGroup )
  n = ncols(G.Z)
  S, a = PolynomialRing(base_field(G), "a" => (1:n, 1:n))
  A = matrix(a)
  return LinearRepresentation( G, A )
end

@Markdown.doc """
    mapping_matrix( rep::LinearRepresentation )

Returns the r×r-matrix A with polynomial entries in the 
variables zᵢⱼ, 1≤i,j≤n of the group G.
"""
function mapping_matrix( rep::LinearRepresentation )
  return rep.A
end

@Markdown.doc """
    group( rep::LinearRepresentation )

Returns the group G which is represented by `rep`.
"""
function group( rep::LinearRepresentation )
  return rep.G
end

@Markdown.doc """
    representation_on_sym_power( rep::LinearRepresentation, d::Int )

Constructs the induced representation on the d-th symmetric power 
Symᵈ V* of the vector space V on which `rep` is defined.
"""
function representation_on_sym_power( rep::LinearRepresentation, d::Int )
  n = ncols(mapping_matrix(rep))	        # the rank of the original representation
  N = binomial( n+d-1, d )	# the rank of the induced representation
  G = group(rep)
  # A polynomial ring containing variables for the coordinates of 
  # the original representation:
  R_z = base_ring(mapping_matrix(rep))
  S, y = PolynomialRing( R_z, [ "y$k" for k in (1:n) ])
  # An empty matrix for the induced representation
  M = MatrixSpace( R_z, N, N )
  B = zero(M)
  # The matrix of the original representation, but promoted to the 
  # new ring S
  @show mapping_matrix(rep)
  @show S
  A = matrix(S.(mapping_matrix(rep)))
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
  return LinearRepresentation( G, B )
end
