@attributes mutable struct SimpleLieAlgebra{C<:RingElement} <: LieAlgebra{C}
#Construct a simple Lie algebra over a given ring with a given root system
  mat_space::MatSpace
  base_ring::Ring
  root_system::RootSystem
  dim::Int
  s::Vector{Symbol}
  root_type::Tuple{String,Int64}

  function SimpleLieAlgebra{C}(R::Ring, S::String, cached::Bool) where {C<:RingElement}
    RS = RootSystem(S)
    di = number_of_roots(RS) + length(RS.simple_roots)
    M = MatrixSpace(R, 1, di)
    s = [Symbol("e_$i") for i in 1:di]
    st=RS.root_system_type
    return get_cached!(SimpleLieAlgebraDict, (R, RS, M), cached) do
    new{C}(M,R,RS,di,s,st)
    end::SimpleLieAlgebra{C}
  end
end

const SimpleLieAlgebraDict = CacheDictType{
  Tuple{Ring,RootSystem,MatSpace},SimpleLieAlgebra
}()

mutable struct SimpleLieAlgebraElem{C<:RingElement} <: LieAlgebraElem{C}
  parent::SimpleLieAlgebra{C}
  mat::MatElem{C}
  function SimpleLieAlgebraElem{C}(
  L::SimpleLieAlgebra{C},
  M::MatElem{C}) where {C<:RingElement}
     return new{C}(L,M)
  end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{SimpleLieAlgebraElem{C}}) where {C<:RingElement} = SimpleLieAlgebra{C}

elem_type(::Type{SimpleLieAlgebra{C}}) where {C<:RingElement} = SimpleLieAlgebraElem{C}

parent(x::SimpleLieAlgebraElem{C}) where {C<:RingElement} = x.parent

base_ring(x::SimpleLieAlgebraElem{C}) where {C<:RingElement} = x.parent.base_ring

base_ring(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.base_ring::parent_type(C)

dim(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.dim

root_system(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.root_system

root_type(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.root_type

matrix(f::SimpleLieAlgebraElem) = f.mat

@doc raw"""
    characteristic(L::SimpleLieAlgebra{T}) -> Int64

Return the characteristic of the base ring of the Lie algebra `L`.
"""
characteristic(L::SimpleLieAlgebra{T}) where T <: RingElement = characteristic(base_ring(L))

# linear indexing 
getindex(f::SimpleLieAlgebraElem, r::Int, c::Int) = f.mat[r, c]

Base.@propagate_inbounds function Base.setindex!(f::SimpleLieAlgebraElem{T}, d, r) where {T <: RingElem}
	matrix(f)[1, r]=d
	return f
end

function symbols(L::SimpleLieAlgebra) 
  return L.s
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::SimpleLieAlgebra)
  print(io, "Simple Lie Algebra over ")
  print(IOContext(io, :compact => true), base_ring(V))
end

show(io::IO, x::SimpleLieAlgebraElem) = show(io, matrix(x))

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Unary operations

function -(f::SimpleLieAlgebraElem)
  R = parent(f)
  return R(-matrix(f))
end

# Binary operations

function +(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  R = parent(f)
  return R(matrix(f) + matrix(g))
end

function -(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  R = parent(f)
  return R(matrix(f) - matrix(g))
end

*(x::RingElem, f::SimpleLieAlgebraElem{T}) where T <: RingElement = parent(f)(x*f.mat)

function *(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  L = parent(f)
  ad = adjoint_matrix(L)
  r = L()
  for i = 1:ncols(matrix(f))
    for j = 1:ncols(matrix(g))
      r = r+f[i]*g[j]*L(transpose(ad[i][:, j]))
   	end
  end
  return r
end

function bracket(x::SimpleLieAlgebraElem{C}, y::SimpleLieAlgebraElem{C}) where {C<:RingElement}
  check_parent(x, y)
  return x*y
end

###############################################################################
#
#   Comparison
#
###############################################################################

function ==(L::SimpleLieAlgebra{T}, M::SimpleLieAlgebra{T}) where T <: RingElement
  (L.root_system).root_system_type != (M.root_system).root_system_type && error("Incompatible Root systems")
  return L.mat_space == M.mat_space
end

function ==(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  return matrix(f) == matrix(g)
end
###############################################################################
#
#   Constructors
#
###############################################################################

@doc raw"""
    zero(L::SimpleLieAlgebra{T}) -> SimpleLieAlgebraElem{T}

Return the zero element of the Lie algebra `L`.
"""
function zero(L::SimpleLieAlgebra{T}) where {T <: RingElement} 
  x = zero(L.mat_space)
	f = SimpleLieAlgebraElem{T}(L, x)
	return f
end

function (L::SimpleLieAlgebra{T})() where {T <: RingElement}
  return zero(L)
end

function (L::SimpleLieAlgebra{T})(A::MatElem{T}) where {T <: RingElem}
  f = SimpleLieAlgebraElem{T}(L, L.mat_space(A))
  return f
end

@doc raw"""
    lie_algebra(R::Ring, S::String, cached::Bool=true) -> SimpleLieAlgebra{elem_type(R)}

Construct the simple Lie algebra over the ring `R` with root system of type `S`
"""
function lie_algebra(R::Ring, S::String, cached::Bool=true)
  T = elem_type(R)
  return SimpleLieAlgebra{T}(R, S, cached)
end

###############################################################################
#
#   Chevalley basis and adjoint matrix
#
###############################################################################

@doc raw"""
    chevalley_basis(L::SimpleLieAlgebra{T}) -> Vector{Vector{SimpleLieAlgebraElem{T}}}

Give the Chevalley basis of the simple Lie algebra `L`
"""
function chevalley_basis(L::SimpleLieAlgebra)
	RS = root_system(L)
	Rp = RS.positive_roots
	sR = RS.simple_roots
	
	n = length(Rp)
	m = length(sR)
	#root vectors
	e1 = elem_type(L)[]
	e2 = elem_type(L)[]
	for i = 1:n
		z1 = L()
		z1[i] = 1
		z2 = L()
		z2[i+n] = 1
		e1 = reduce(vcat, (e1, z1))
		e2 = reduce(vcat, (e2, z2))
	end
	#basis for cartan algebra
	e3 = elem_type(L)[]
	for i = 1:m
		z = L()
		z[i+2*n] = 1
		e3 = reduce(vcat, (e3, z))
	end
	return [e1, e2, e3]
end 

@doc raw"""
    adjoint_matrix(L::SimpleLieAlgebra{T}) -> Vector{Any}

Give the adjoint matrices of all basis vectors acting on the Lie algebra `L` with respect 
to the Chevalley basis of `L`
"""
function adjoint_matrix(L::SimpleLieAlgebra{T}) where T <: RingElement #computes the adjoint matrix with respect to the Chevalley basis.
	St = L.root_type
	R = base_ring(L)	
	n = St[2]
	
	Q = GAP.Globals.Rationals
	S = GAP.Obj(St[1])
	LG = GAP.Globals.SimpleLieAlgebra(S, n, Q) #define the Lie algebra in Gap over the rationals, we are interested in the structure constants
  #they remain the same over rings of characteristic 0, otherwise we have the structure constants mod char(R)
	d = GAP.Globals.Dimension(LG)
	
	ch = GAP.Globals.ChevalleyBasis(LG)
	ch1 = GAP.Globals.ShallowCopy(ch[1])
	ch2 = GAP.Globals.ShallowCopy(ch[2])
	ch3 = GAP.Globals.ShallowCopy(ch[3])
	ch = ch1
	GAP.Globals.Append(ch, ch2)
	GAP.Globals.Append(ch, ch3)
	BL = GAP.Globals.Basis(LG, ch)
	ad = dense_matrix_type(R)[]
	for i = 1:d
		Ad = GAP.Globals.AdjointMatrix(BL, BL[i])
		A = [[Ad[i][j] for j = 1:length(Ad[i])] for i = 1:length(Ad)]
		MA = zero_matrix(R, d, d)
		for j = 1:d
			MA[j, 1:d] = A[j]
		end
		push!(ad, MA)
		
	end 
	return ad
end

