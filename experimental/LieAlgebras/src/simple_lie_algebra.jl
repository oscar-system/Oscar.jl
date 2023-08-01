@attributes mutable struct SimpleLieAlgebra{C<:RingElement} <: LieAlgebra{C}
#Construct a simple Lie algebra over a given ring with a given root system
  base_ring::Ring
  root_system::RootSystem
  dim::Int
  s::Vector{Symbol}
  root_type::Tuple{Symbol,Int64}
  struct_consts::Matrix{SRow{C}}
  function SimpleLieAlgebra{C}(
    R::Ring, S::Symbol, n::Int64; cached::Bool=true
  ) where {C<:RingElement}
    RS = root_system(S,n)
    dimL = number_of_roots(RS) + length(RS.simple_roots)
    s = [Symbol("e_$i") for i in 1:dimL]
    st = root_system_type(RS)
    #get the structure constants of the Lie algebra L 
    #note that it is enoough to do this over Q, as we can later coerce the constants
    #into the ring R
    coeffs_iso = inv(Oscar.iso_oscar_gap(QQ))
    LG = GAP.Globals.SimpleLieAlgebra(
      GAP.Obj(S), n, domain(coeffs_iso)
    )
    sc_table_G =
      (
        entry -> (entry[1], Vector{elem_type(RO)}(map(coeffs_iso, entry[2])))
      ).(
        Matrix{Tuple{Vector{Int},Vector{GAP.Obj}}}(
          (GAP.Globals.StructureConstantsTable(GAPWrap.Basis(LG)))[1:dimL]
        )
      )
    struct_consts = Matrix{SRow{elem_type(RO)}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
      struct_consts[i, j] = sparse_row(
        RO, Tuple{Int,elem_type(R)}[(k, R(c)) for (k, c) in zip(sc_table_G[i, j]...)]
      )
    end
    return get_cached!(SimpleLieAlgebraDict, (R, RS), cached) do
    new{C}(R, RS, dimL, s, st, struct_consts)
    end::SimpleLieAlgebra{C}
  end
end

const SimpleLieAlgebraDict = CacheDictType{
  Tuple{Ring,RootSystem},SimpleLieAlgebra
}()

mutable struct SimpleLieAlgebraElem{C<:RingElement} <: LieAlgebraElem{C}
  parent::SimpleLieAlgebra{C}
  mat::MatElem{C}
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{SimpleLieAlgebraElem{C}}) where {C<:RingElement} = SimpleLieAlgebra{C}

elem_type(::Type{SimpleLieAlgebra{C}}) where {C<:RingElement} = SimpleLieAlgebraElem{C}

parent(x::SimpleLieAlgebraElem{C}) where {C<:RingElement} = x.parent

base_ring(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.base_ring::parent_type(C)

dim(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.dim

root_system(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.root_system

root_type(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.root_type

#  Base.@propagate_inbounds function Base.setindex!(f::SimpleLieAlgebraElem{T}, d, r) where {T <: RingElem}
#    Generic._matrix(f)[1, r] = d
#    return f
#  end

function symbols(L::SimpleLieAlgebra) 
  return L.s
end

###############################################################################
#
#   String I/O
#
###############################################################################
function Base.show(io::IO, ::MIME"text/plain", L::SimpleLieAlgebra)
  io = pretty(io)
  println(io, "Simple Lie algebra")
  println(
    io, Indent(), "of type $(string(root_type(L)[1])*string(root_type(L)[2]))", Dedent()
  )
  println(io, Indent(), "of dimension $(dim(L))", Dedent())
  print(io, "over ")
  print(io, Lowercase(), coefficient_ring(L))
end

function Base.show(io::IO, L::SimpleLieAlgebra)
  if get(io, :supercompact, false)
    print(io, "Simple Lie algebra")
  else
    io = pretty(io)
    print(io, "Simple Lie algebra over ", Lowercase())
    print(IOContext(io, :supercompact => true), coefficient_ring(L))
  end
end

# function Base.show(io::IO, V::SimpleLieAlgebra)
#   print(io, "Simple Lie Algebra over ")
#   print(IOContext(io, :compact => true), base_ring(V))
# end


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Binary operations

function bracket(
  x::SimpleLieAlgebraElem{C},
  y::SimpleLieAlgebraElem{C},
) where {C<:RingElement}
  check_parent(x, y)
  L = parent(x)
  mat = sum(
    cxi * cyj * L.struct_consts[i, j] for (i, cxi) in enumerate(coefficients(x)),
    (j, cyj) in enumerate(coefficients(y))
  )
  return L(mat)
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(L::SimpleLieAlgebra{T}, M::SimpleLieAlgebra{T}) where T <: RingElement
  return base_ring(L) == base_ring(M) &&
         dim(L) == dim(M) &&
         root_system(L) == root_system(M)
end

function Base.hash(x::SimpleLieAlgebra{C}, h::UInt) where {C<:RingElement}
  b = 0x9072419641dbc7d0 % UInt
  h = hash(base_ring(L), h)
  h = hash(dim(L), h)
  h = hash(mat_space(L),h)
  h = hash(root_system(L),h)
  h = hash(root_type(L),h)
  return xor(h, b)
end

###############################################################################
#
#   Constructors
#
###############################################################################

@doc raw"""
    lie_algebra(R::Ring, S::String, cached::Bool=true) -> SimpleLieAlgebra{elem_type(R)}

Construct the simple Lie algebra over the ring `R` with root system of type `S`
The internally used basis of this Lie algebra is the Chevalley basis.
"""
function lie_algebra(R::Ring, S::Symbol, n::Int64; cached::Bool=true)
  @req S in [:A, :B, :C, :D, :E, :F, :G] "Unknown Dynkin type"
  T = elem_type(R)
  return SimpleLieAlgebra{T}(R, S, n; cached)
end

###############################################################################
#
#   Chevalley basis and adjoint matrix
#
###############################################################################

@doc raw"""
    chevalley_basis(L::SimpleLieAlgebra{T}) -> Vector{Vector{SimpleLieAlgebraElem{T}}}

Give the Chevalley basis of the simple Lie algebra `L` in three vectors, stating first the positive root vectors, 
then the negative root vectors and finally the basis of the Cartan subalgebra. The order of root vectors corresponds
to the order or the roots in the root system.
"""
function chevalley_basis(L::SimpleLieAlgebra)
  RS = root_system(L)
  Rp = RS.positive_roots
  sR = RS.simple_roots
  n = length(Rp)
  m = length(sR)
  B = basis(L)
  #root vectors
  r_plus = [B[i] for i = 1:n]
  r_minus = [B[i] for i = (n + 1):(2 * n)]
  #basis for cartan algebra
  h = [B[i] for i = (2 * n + 1): dim(L)]
  return [r_plus, r_minus, h]
end 

@doc raw"""
    adjoint_matrix(x::SimpleLieAlgebraElem{T}) -> MatSpaceElem{T}

Give the adjoint matrix of the Lie algebra element `x` with respect 
to the Chevalley basis.
The adjoint matrix is the matrix of the adjoint representation of the element x 
w.r.t. the Chevalley basis B. The adjoint map is the left multiplication by x. 
The j-th column of the resulting matrix represents the image of the the j-th basis vector of B 
under left multiplication by x. 
"""
function adjoint_matrix(x::SimpleLieAlgebraElem{T}) where {T<:RingElement} 
  #computes the adjoint matrix of the element x with respect to the chevalley basis.
  L = parent(x)
  #compute x acting on the basis vectors
  Ad = zero_matrix(base_ring(x), dim(L), dim(L))
  for i in 1:dim(L)
    Ad[:,i]=coefficients(bracket(x,basis(L,i)))
  end
  return Ad
end

@doc raw"""
    adjoint_matrix(x::SimpleLieAlgebraElem{T}) -> MatSpaceElem{T}

Give the adjoint matrix of the Lie algebra element `x` with respect 
to the basis `B`. Note that we expect the elements in `B` to be given with coefficients w.r.t 
the Chevalley basis but the coefficents of `x` are given w.r.t B.
The adjoint matrix is the matrix of the adjoint representation of the element x 
w.r.t. the `B`. The adjoint map is the left multiplication by x. 
The j-th column of the resulting matrix represents the image of the the j-th basis vector of B 
under left multiplication by x. 
"""
function adjoint_matrix(
  x::SimpleLieAlgebraElem{T}, B::Vector{SimpleLieAlgebraElem{T}}
) where {T<:RingElement}
  R = base_ring(x)
  #compute the base change matrix from the Chevalley basis to B 
  M_B = matrix(R, [[B[i][j] for j in 1:length(Generic._matrix(B[i]))] for i in 1:length(B)])
  iM_B = inv(M_B)
  #transform x to coordinates over Chevalley basis
  y = L(Generic._matrix(x) * M_B)
  ad_ch = adjoint_matrix(y)
  ad = iM_B * ad_ch * M_B
  return ad
end

@doc raw"""
    adjoint_matrix(L::SimpleLieAlgebra{T}) -> Vector{MatSpaceElem{T}}

Give the adjoint matrices of all basis vectors acting on the Lie algebra `L` with respect 
to the Chevalley basis of `L`.
The i-th adjoint matrix is the matrix of the adjoint representation of the i-th basis element x 
w.r.t. the Chevalley basis B. The adjoint map is the left multiplication by x. 
The j-th column of the resulting matrix represents the image of the the j-th basis vector of B 
under left multiplication by x. 
"""
function adjoint_matrix(L::SimpleLieAlgebra{T}) where T <: RingElement #computes the adjoint matrix with respect to the Chevalley basis.
  R = base_ring(L)
  ad = dense_matrix_type(R)[]
  d = dim(L)
  B = basis(L)
  for i = 1:d
    Ad = adjoint_matrix(B[i])
    push!(ad, Ad)
  end 
  return ad
end
