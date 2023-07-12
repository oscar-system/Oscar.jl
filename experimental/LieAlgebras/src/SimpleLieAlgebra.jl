GAP.Packages.load("sla")
using Oscar
using .Generic

import Oscar: parent_type, elem_type, base_ring, parent, MatElem, symbols, dim

import .Generic: CacheDictType

import AbstractAlgebra: get_cached!

import Base: show, +, -, *, ^, ==, !=, inv, isone, iszero, one, zero, rand, deepcopy_internal, hash, setindex!


@attributes mutable struct SimpleLieAlgebra{C<:RingElement} <: LieAlgebra{C}
#Construct a simple Lie algebra over a given ring with a given root system
  mat_space::MatSpace
  base_ring::Ring
  root_system::RootSystem
  dim::Int
  s::Vector{Symbol}
  root_type::Vector{Union{String, Int64}}

  function SimpleLieAlgebra{C}(
    R::Ring,
    S::String,
    cached::Bool
  ) where {C<:RingElement}
    l = length(S)
    n = parse(Int64, S[2:l])
    Q = GAP.Globals.Rationals
    S1 = GAP.Obj(string(S[1]))
    sL = GAP.Globals.SimpleLieAlgebra(S1, n, Q)
    di = GAP.Globals.Dimension(sL)
    M = MatrixSpace(R, 1, di)
    RS = RootSystem(S)
    s = [Symbol("e_$i") for i in 1:di]
    rt=Union{String, Int64}[S[1:1], n]
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
  base_ring::Ring
  function SimpleLieAlgebraElem{C}(
  L::SimpleLieAlgebra{C},
  M::MatElem{C}) where {C<:RingElement}
     return new{C}(L,M,base_ring(L))
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

base_ring(x::SimpleLieAlgebraElem{C}) where {C<:RingElement} = x.base_ring

base_ring(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.base_ring::parent_type(C)

dim(L::SimpleLieAlgebra{C}) where {C<:RingElement} = L.dim

characteristic(R::SimpleLieAlgebra{T}) where T <: RingElement = characteristic(base_ring(R))

# linear indexing 
getindex(f::SimpleLieAlgebraElem, r::Int, c::Int) = f.mat[r, c]

Base.@propagate_inbounds function Base.setindex!(f::SimpleLieAlgebraElem{T}, d, r) where {T <: RingElem}
	f.mat[1, r]=d
	return f
end

###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::SimpleLieAlgebra{C}) where {C<:RingElement}
  print(io, "Simple Lie Algebra over ")
  print(IOContext(io, :compact => true), base_ring(V))
end

show(io::IO, x::SimpleLieAlgebraElem{C}) where {C<:RingElement} = show(io, x.mat)

function symbols(L::SimpleLieAlgebra{C}) where {C<:RingElement}
  return L.s
end

###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Unary operations

function -(f::SimpleLieAlgebraElem)
  R = parent(f)
  return R(-f.mat)
end

# Binary operations

function +(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  R = parent(f)
  return R(f.mat + g.mat)
end

function -(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  R = parent(f)
  return R(f.mat - g.mat)
end

*(x::RingElem, f::SimpleLieAlgebraElem{T}) where T <: RingElement = parent(f)(x*f.mat)

function *(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  L = parent(f)
  ad = AdjointMatrix(L)
  r = L()
  for i = 1:ncols(f.mat)
    for j = 1:ncols(g.mat)
      r = r+f[i]*g[j]*L(transpose(ad[i][:, j]))
   	end
  end
  return r
end

function bracket(
  x::SimpleLieAlgebraElem{C}, y::SimpleLieAlgebraElem{C}
) where {C<:RingElement}
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

function !=(L::SimpleLieAlgebra{T}, M::SimpleLieAlgebra{T}) where T <: RingElement
  if (L.root_system).root_system_type != (M.root_system).root_system_type 
    return true
  else
   	return L.mat_space != M.mat_space
  end
end

function ==(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  return f.mat == g.mat
end

function !=(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  return f.mat != g.mat
end

function isequal(f::SimpleLieAlgebraElem{T}, g::SimpleLieAlgebraElem{T}) where T <: RingElement
  parent(f) != parent(g) && error("Incompatible Lie algebras")
  return isequal(f.mat, g.mat)
end

function in(f::SimpleLieAlgebraElem, L::SimpleLieAlgebra)
	if f.parent == L
    return true
	else
		return false
	end
end
###############################################################################
#
#   Constructors
#
###############################################################################

function zero(R::SimpleLieAlgebra{T}) where {T <: RingElement} 
  x = zero(R.mat_space)
	f = SimpleLieAlgebraElem{T}(R, x)
	return f
end

function (R::SimpleLieAlgebra{T})() where {T <: RingElement}
  return zero(R)
end

function (R::SimpleLieAlgebra{T})(A::MatElem{T}) where {T <: RingElem}
  C = elem_type(R)
  f = SimpleLieAlgebraElem{T}(R, R.mat_space(A))
  return f
end

@doc raw"""
    LieAlgebra(R::Ring, S::String, cached::Bool=true) -> SimpleLieAlgebra{elem_type(R)}
Construct the simple Lie algebra over the ring `R` with root system of type `S`
"""
function LieAlgebra(R::Ring, S::String, cached::Bool=true)
  T = elem_type(R)
  return SimpleLieAlgebra{T}(R, S, cached)
end

###############################################################################
#
#   Chevalley basis and adjoint matrix
#
###############################################################################

@doc raw"""
    ChevalleyBasis(L::SimpleLieAlgebra{T}) -> Vector{Vector{SimpleLieAlgebraElem{T}}}
Give the Chevalley basis of a simple Lie algebra
"""
function ChevalleyBasis(L::SimpleLieAlgebra{T}) where {T<: RingElement}
	RS = L.root_system
	Rp = RS.positive_roots
	sR = RS.simple_roots
	
	n = length(Rp)
	m = length(sR)
	#root vectors
	e1 = [ ]
	e2 = [ ]
	for i = 1:n
		z1 = L()
		z1[i] = 1
		z2 = L()
		z2[i+n] = 1
		e1 = reduce(vcat, (e1, z1))
		e2 = reduce(vcat, (e2, z2))
	end
	#basis for cartan algebra
	e3 = [ ]
	for i = 1:m
		z = L()
		z[i+2*n] = 1
		e3 = reduce(vcat, (e3, z))
	end
	return [e1, e2, e3]
end 

@doc raw"""
    AdjointMatrix(L::SimpleLieAlgebra{T}) -> Vector{Any}
Give the Adjoint matrices of all basis vectors acting on the Lie algebra L with respect to the Chevalley basis of L
"""
function AdjointMatrix(L::SimpleLieAlgebra{T}) where T <: RingElement #computes the adjoint matrix with respect to the Chevalley basis.
	St = L.root_type
	R = base_ring(L)	
	n = St[2]
	
	Q = GAP.Globals.Rationals
	S = GAP.Obj(St[1])
	LG = GAP.Globals.SimpleLieAlgebra(S, n, Q) #define the Lie algebra in Gap over the rationals
	d = GAP.Globals.Dimension(LG)
	
	ch = GAP.Globals.ChevalleyBasis(LG)
	ch1 = GAP.Globals.ShallowCopy(ch[1])
	ch2 = GAP.Globals.ShallowCopy(ch[2])
	ch3 = GAP.Globals.ShallowCopy(ch[3])
	ch = ch1
	GAP.Globals.Append(ch, ch2)
	GAP.Globals.Append(ch, ch3)
	BL = GAP.Globals.Basis(LG, ch)
	M = MatrixSpace(R, d, d)
	ad = [ ]
	for i = 1:d
		Ad = GAP.Globals.AdjointMatrix(BL, BL[i])
		A = [[Ad[i][j] for j = 1:length(Ad[i])] for i = 1:length(Ad)]
		MA = M()
		for j = 1:d
			MA[j, 1:d] = A[j]
		end
		append!(ad,[MA])
		
	end 
	return ad
end

