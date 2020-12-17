import AbstractAlgebra: MatElem, matrix, MatSpace, parent_type, Ring, RingElem
import Hecke: base_ring, det, fq_nmod, FqNmodFiniteField, nrows, tr, trace
import GAP: FFE

export
    general_linear_group,
    matrix_group,
    MatrixGroup,
    MatrixGroupElem,
    omega_group,
    orthogonal_group,
    special_linear_group,
    special_orthogonal_group,
    special_unitary_group,
    symplectic_group,
    unitary_group,
    GL, GO, GU, SL, SO, Sp, SU



abstract type AbstractMatrixGroupElem <: GAPGroupElem{GAPGroup} end

# NOTE: always defined are deg, ring and at least one between { X, gens, descr }
# NOTE: the fields ring_iso and mat_iso are always defined if the field X is
"""
    MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup
Type of groups `G` of `n x n` matrices over the ring `R`, where `n = degree(G)` and `R = base_ring(G)`.

At the moment, only rings of type `FqNmodFiniteField` are supported.
"""
mutable struct MatrixGroup{RE<:RingElem, T<:MatElem{RE}} <: GAPGroup
   deg::Int
   ring::Ring
   X::GapObj
   gens::Vector{<:AbstractMatrixGroupElem}
   descr::Symbol                       # e.g. GL, SL, symbols for isometry groups
   ring_iso::GenRingIso
   mat_iso::GenMatIso
   AbstractAlgebra.@declare_other

   MatrixGroup(m::Int, F::Ring) = new{elem_type(F), dense_matrix_type(elem_type(F))}(m,F)

end

# build a MatrixGroup given the underlying GAP object
function MatrixGroup(m::Int, F::S, G_gap::GapObj) where S<:Ring
   r = MatrixGroup(m,F)
   r.X = G_gap
   return r
end

# build a MatrixGroup given a list of generators
function MatrixGroup(m::Int, F::Ring, L::Vector)
   r = MatrixGroup(m,F)
   r.gens = L
   for l in r.gens l.parent = r end
   return r
end


# NOTE: at least one of the fields :elm and :X must always defined, but not necessarily both of them.
"""
    MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}} <: GAPGroupElem{MatrixGroup}

Elements of a group of type `MatrixGroup{RE<:RingElem, T<:MatElem{RE}}`
"""
mutable struct MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}} <: AbstractMatrixGroupElem
   parent::MatrixGroup{RE, T}
   elm::T                         # Oscar matrix
   X::GapObj                     # GAP matrix. If x isa MatrixGroupElem, then x.X = x.parent.mat_iso(x.elm)

   #MatrixGroupElem(G::MatrixGroup{RE,MatElem{RE}}, x::MatElem{RE}) where RE <: RingElem = new{RE, MatElem{RE}}(G,x)
   MatrixGroupElem(G::MatrixGroup, x::MatElem) = new{elem_type(G.ring),dense_matrix_type(elem_type(G.ring))}(G,x)

   #function MatrixGroupElem(G::MatrixGroup{RE,MatElem{RE}}, x_gap::GapObj) where RE <: RingElem
   function MatrixGroupElem(G::MatrixGroup, x_gap::GapObj)
      z = new{elem_type(G.ring),dense_matrix_type(elem_type(G.ring))}()
      z.parent = G
      z.X = x_gap
      return z
   end

end

# build a MatrixGroupElem given its underlying GAP object
# WARNING: this does not check whether the element actually lies in the group G
function MatrixGroupElem(G::MatrixGroup, x::MatElem, x_gap::GapObj)
   z = MatrixGroupElem(G,x)
   z.X = x_gap
   return z
end


ring_elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = S
mat_elem_type(::Type{MatrixGroup{S,T}}) where {S,T} = T
_gap_filter(::Type{<:MatrixGroup}) = GAP.Globals.IsMatrixGroup

########################################################################
#
# Basic
#
########################################################################

function Base.show(io::IO, x::MatrixGroup)
   AbstractAlgebra.@show_name(io, x)
   AbstractAlgebra.@show_special(io,x)

   if isdefined(x, :descr)
      if x.descr==:GU || x.descr==:SU
         print(io, string(x.descr), "(",x.deg,",",characteristic(x.ring)^(div(degree(x.ring),2)),")")
      else
         print(io, string(x.descr), "(",x.deg,",",order(x.ring),")")
      end
   else
      print(io, "Matrix group of degree ", x.deg, " over ")
      show(IOContext(io, :compact => true), x.ring)
   end
end

function Base.show(io::IO, x::MatrixGroupElem)          #TODO: is this correct?
   if isdefined(x, :elm)
      show(io, "text/plain", x.elm)
#      print(io, x.elm)
   else
      print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.X)))
   end
end

group_element(G::MatrixGroup, x::GapObj) = MatrixGroupElem(G,x)

function assign_from_description(G::MatrixGroup)
   if G.descr==:GL G.X=GAP.Globals.GL(G.deg,G.ring_iso.codomain)
   elseif G.descr==:SL G.X=GAP.Globals.SL(G.deg,G.ring_iso.codomain)
   elseif G.descr==:Sp G.X=GAP.Globals.Sp(G.deg,G.ring_iso.codomain)
   elseif G.descr==Symbol("GO+") G.X=GAP.Globals.GO(1,G.deg,G.ring_iso.codomain)
   elseif G.descr==Symbol("SO+") G.X=GAP.Globals.SO(1,G.deg,G.ring_iso.codomain)
   elseif G.descr==Symbol("Omega+") G.X=GAP.Globals.Omega(1,G.deg,Int(order(G.ring)))
   elseif G.descr==Symbol("GO-") G.X=GAP.Globals.GO(-1,G.deg,G.ring_iso.codomain)
   elseif G.descr==Symbol("SO-") G.X=GAP.Globals.SO(-1,G.deg,G.ring_iso.codomain)
   elseif G.descr==Symbol("Omega-") G.X=GAP.Globals.Omega(-1,G.deg,Int(order(G.ring)))
   elseif G.descr==:GO G.X=GAP.Globals.GO(0,G.deg,G.ring_iso.codomain)
   elseif G.descr==:SO G.X=GAP.Globals.SO(0,G.deg,G.ring_iso.codomain)
   elseif G.descr==:Omega G.X=GAP.Globals.Omega(0,G.deg,Int(order(G.ring)))
   elseif G.descr==:GU G.X=GAP.Globals.GU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) ))
   elseif G.descr==:SU G.X=GAP.Globals.SU(G.deg,Int(characteristic(G.ring)^(div(degree(G.ring),2) ) ))
   else error("unsupported description")
   end
end

# return the G.sym if isdefined(G, :sym); otherwise, the field :sym is computed and set using information from other defined fields
# NOTE: if G.X has to be set, then also the fields G.ring_iso and G.mat_iso are set
function Base.getproperty(G::MatrixGroup, sym::Symbol)

   if isdefined(G,sym) return getfield(G,sym) end

   if sym === :mat_iso || sym === :ring_iso
      fm = gen_mat_iso(G.deg, G.ring)
      G.ring_iso=fm.fr
      G.mat_iso=fm

   elseif sym === :X
      if isdefined(G,:descr)
         if !isdefined(G,:ring_iso)
            fm = gen_mat_iso(G.deg, G.ring)
            G.ring_iso=fm.fr
            G.mat_iso=fm
         end
         assign_from_description(G)
      elseif isdefined(G,:gens)
         if !isdefined(G,:ring_iso)
            fm = gen_mat_iso(G.deg, G.ring)
            G.ring_iso=fm.fr
            G.mat_iso=fm
         end
         V = GAP.julia_to_gap([g.X for g in G.gens])
         G.X=GAP.Globals.Group(V)
      end

   elseif sym === :gens
      if !isdefined(G, :X)
         if !isdefined(G,:ring_iso)
            fm = gen_mat_iso(G.deg, G.ring)
            G.ring_iso=fm.fr
            G.mat_iso=fm
         end
         assign_from_description(G)
      end
      L = GAP.Globals.GeneratorsOfGroup(getfield(G,:X))
      G.gens=[MatrixGroupElem(G,G.mat_iso(L[i]),L[i]) for i in 1:length(L)]

   end

   return getfield(G, sym)

end


function Base.getproperty(x::MatrixGroupElem, sym::Symbol)

   if isdefined(x,sym) return getfield(x,sym) end
   if sym === :X && !isdefined(x,:X)
      if !isdefined(x.parent,:ring_iso)
         fm = gen_mat_iso(x.parent.deg, x.parent.ring)
         setfield!(x.parent,:ring_iso,fm.fr)
         setfield!(x.parent,:mat_iso,fm)
      end
      setfield!(x, :X, x.parent.mat_iso(x.elm))
   elseif sym == :elm && !isdefined(x, :elm)                        # assuming that x.X, hence x.parent.ring_iso, are defined
      setfield!(x, :elm, x.parent.mat_iso(x.X))
   end
   return getfield(x,sym)
end

Base.IteratorSize(::Type{<:MatrixGroup}) = Base.SizeUnknown()

function Base.iterate(G::MatrixGroup)
  L=GAP.Globals.Iterator(G.X)
  @assert ! GAP.Globals.IsDoneIterator(L)
  i = GAP.Globals.NextIterator(L)
  return MatrixGroupElem(G, G.mat_iso(i),i), L
end

function Base.iterate(G::MatrixGroup, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return MatrixGroupElem(G, G.mat_iso(i),i), state
end

########################################################################
#
# Membership
#
########################################################################


function ==(G::MatrixGroup,H::MatrixGroup)
   G.deg==H.deg || return false
   G.ring==H.ring || return false
   if isdefined(G, :descr) && isdefined(H, :descr)
      return G.descr == H.descr
   end
   if isdefined(G, :gens) && isdefined(H, :gens)
      if G.gens==H.gens return true  end       #TODO this works only if the == between elements does not keep track of the parent
   end
   return G.X==H.X
end


# this saves the value of x.X
# x_gap = x.X if this is already known, x_gap = nothing otherwise
function lies_in(x::MatElem, G::MatrixGroup, x_gap)
   if base_ring(x)!=G.ring || nrows(x)!=G.deg return false, x_gap end
   if isone(x) return true, x_gap end
   if isdefined(G,:gens)
      for g in G.gens
         if x==g.elm
            return true, x_gap
         end
      end
   end
   if isdefined(G,:descr) && G.descr==:GL
      if det(x)!=0
         return true, x_gap
      else
         return false, x_gap
      end
   elseif isdefined(G,:descr) && G.descr==:SL
      if det(x)==1
         return true, x_gap
      else
         return false, x_gap
      end
   else
      if x_gap==nothing x_gap = G.mat_iso(x) end
     # x_gap !=nothing || x_gap = G.mat_iso(x)
      return GAP.Globals.in(x_gap,G.X), x_gap
   end
end

Base.in(x::MatElem, G::MatrixGroup) = lies_in(x,G,nothing)[1]

function Base.in(x::MatrixGroupElem, G::MatrixGroup)
   if isdefined(x,:X) return lies_in(x.elm,G,x.X)[1]
   else return lies_in(x.elm,G,nothing)[1]
   end
end

# embedding an element of type MatElem into a group G
# if check=false, there are no checks on the condition `x in G`
function (G::MatrixGroup)(x::MatElem; check=true)
   if check
      vero, x_gap = lies_in(x,G,nothing)
      vero || throw(ArgumentError("Element not in the group"))
      if x_gap==nothing return MatrixGroupElem(G,x)
      else return MatrixGroupElem(G,x,x_gap)
      end
   else
      return MatrixGroupElem(G,x)
   end
end

# embedding an element of type MatrixGroupElem into a group G
# if check=false, there are no checks on the condition `x in G`
function (G::MatrixGroup)(x::MatrixGroupElem; check=true)
   if !check
      z = x
      z.parent = G
      return z
   end
   if isdefined(x,:X)
      if isdefined(x,:elm)
         vero = lies_in(x.elm,G,x.X)[1]
         vero || throw(ArgumentError("Element not in the group"))
         return MatrixGroupElem(G,x.elm,x.X)
      else
         GAP.Globals.in(x.X, G.X) || throw(ArgumentError("Element not in the group"))
         return MatrixGroupElem(G,x.X)
      end
   else
      vero, x_gap = lies_in(x.elm,G,nothing)
      vero || throw(ArgumentError("Element not in the group"))
      if x_gap==nothing return MatrixGroupElem(G,x.elm)
      else return MatrixGroupElem(G,x.elm,x_gap)
      end
   end
end

# embedding a nxn array into a group G
function (G::MatrixGroup)(L::AbstractVecOrMat; check=true)
   x = matrix(G.ring, G.deg, G.deg, L)
   return G(x; check=check)
end


########################################################################
#
# Methods on elements
#
########################################################################

# we are not currently keeping track of the parent; two elements coincide iff their matrices coincide
function ==(x::MatrixGroupElem{S,T},y::MatrixGroupElem{S,T}) where {S,T}
   if isdefined(x,:X) && isdefined(y,:X) return x.X==y.X
   else return x.elm==y.elm
   end
end

# if the parents are different, the parent of the output product is set as GL(n,q)
function _prod(x::MatrixGroupElem,y::MatrixGroupElem)
   G = x.parent==y.parent ? x.parent : GL(x.parent.deg, x.parent.ring)
   if isdefined(x,:X) && isdefined(y,:X) && !(isdefined(x,:elm) && isdefined(y,:elm))
      return MatrixGroupElem(G, x.X*y.X)
   else
      return MatrixGroupElem(G, x.elm*y.elm)
   end
end

# Base.:* is defined in src/Groups/GAPGroups.jl

Base.:*(x::MatrixGroupElem, y::fq_nmod_mat) = x.elm*y
Base.:*(x::fq_nmod_mat, y::MatrixGroupElem) = x*y.elm

Base.:^(x::MatrixGroupElem, n::Int) = MatrixGroupElem(x.parent, x.elm^n)

Base.isone(x::MatrixGroupElem) = isone(x.elm)

Base.inv(x::MatrixGroupElem) = MatrixGroupElem(x.parent, inv(x.elm))

# if the parents are different, the parent of the output is set as GL(n,q)
function Base.:^(x::MatrixGroupElem, y::MatrixGroupElem)
   G = x.parent==y.parent ? x.parent : GL(x.parent.deg, x.parent.ring)
   if isdefined(x,:X) && isdefined(y,:X) && !(isdefined(x,:elm) && isdefined(y,:elm))
      return MatrixGroupElem(G, inv(y.X)*x.X*y.X)
   else
      return MatrixGroupElem(G,inv(y.elm)*x.elm*y.elm)
   end
end

comm(x::MatrixGroupElem, y::MatrixGroupElem) = inv(x)*conj(x,y)

"""
    det(x::MatrixGroupElem)
    
Return the determinant of `x`.
"""
det(x::MatrixGroupElem) = det(x.elm)

"""
    base_ring(x::MatrixGroupElem)

Return the base ring of `x`.
"""
base_ring(x::MatrixGroupElem) = x.parent.ring

parent(x::MatrixGroupElem) = x.parent

"""
    matrix(x::MatrixGroupElem)

Return the underlying `AbstractAlgebra` matrix of `x`.
"""
matrix(x::MatrixGroupElem) = x.elm

Base.getindex(x::MatrixGroupElem, i::Int, j::Int) = x.elm[i,j]

"""
    nrows(x::MatrixGroupElem)

Return the number of rows of the given matrix.
"""
nrows(x::MatrixGroupElem) = x.parent.deg

"""
    trace(x::MatrixGroupElem)
    tr(x::MatrixGroupElem)

Return the trace of `x`.
"""
trace(x::MatrixGroupElem) = trace(x.elm)
tr(x::MatrixGroupElem) = tr(x.elm)

order(x::MatrixGroupElem) = GAP.Globals.Order(x.X)

#FIXME for the following functions, the output may not belong to the parent group of x
#=
frobenius(x::MatrixGroupElem, n::Int) = MatrixGroupElem(x.parent, matrix(x.parent.ring, x.parent.deg, x.parent.deg, [frobenius(y,n) for y in x.elm]))
frobenius(x::MatrixGroupElem) = frobenius(x,1)

transpose(x::MatrixGroupElem) = MatrixGroupElem(x.parent, transpose(x.elm))
=#

########################################################################
#
# Methods on groups
#
########################################################################

"""
    base_ring(G::MatrixGroup)

Return the base ring of the matrix group `G`.
"""
base_ring(G::MatrixGroup) = G.ring

"""
    degree(G::MatrixGroup)

Return the degree of the matrix group `G`, i.e. the number of rows of its matrices.
"""
degree(G::MatrixGroup) = G.deg

Base.one(G::MatrixGroup) = MatrixGroupElem(G, identity_matrix(G.ring, G.deg))

function Base.rand(G::MatrixGroup)
   x_gap = GAP.Globals.Random(G.X)
   x_oscar = G.mat_iso(x_gap)
   return MatrixGroupElem(G,x_oscar,x_gap)
end

gens(G::MatrixGroup) = G.gens

gen(G::MatrixGroup, i::Int) = gens(G)[i]

Base.getindex(G::MatrixGroup, i::Int) = gen(G, i)

ngens(G::MatrixGroup) = length(G.gens)

order(G::MatrixGroup) = GAP.Globals.Order(G.X)

########################################################################
#
# Constructors
#
########################################################################

"""
    general_linear_group(n::Int, q::Int)
    general_linear_group(n::Int, F::FqNmodFiniteField)
    GL = general_linear_group

Return the general linear group of dimension `n` either over the field `F` or the field `GF(q)`.
"""
function general_linear_group(n::Int, F::Ring)
   G = MatrixGroup(n,F)
   G.descr = :GL
   return G
end

function general_linear_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return general_linear_group(n, GF(b,a)[1])
end

function general_linear_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.GL(_gap_filter(T), n, q))
   elseif T<:MatrixGroup   return general_linear_group(n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

function general_linear_group(::Type{T}, n::Int, F::Ring) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.GL(_gap_filter(T), n, Int(order(F))))
   elseif T<:MatrixGroup   return general_linear_group(n,F)
   else throw(ArgumentError("Type not supported"))
   end
end

"""
    special_linear_group(n::Int, q::Int)
    special_linear_group(n::Int, F::FqNmodFiniteField)
    SL = special_linear_group

Return the special linear group of dimension `n` either over the field `F` or the field `GF(q)`.
"""
function special_linear_group(n::Int, F::Ring)
   G = MatrixGroup(n,F)
   G.descr = :SL
   return G
end

function special_linear_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return special_linear_group(n, GF(b,a)[1])
end

function special_linear_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.SL(_gap_filter(T), n, q))
   elseif T<:MatrixGroup   return special_linear_group(n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

function special_linear_group(::Type{T}, n::Int, F::Ring) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.SL(_gap_filter(T), n, Int(order(F))))
   elseif T<:MatrixGroup   return special_linear_group(n,F)
   else throw(ArgumentError("Type not supported"))
   end
end


"""
    symplectic_group(n::Int, q::Int)
    symplectic_group(n::Int, F::FqNmodFiniteField)
    Sp = symplectic_group

Return the symplectic group of dimension `n` either over the field `F` or the field `GF(q)`. The dimension `n` must be even.
"""
function symplectic_group(n::Int, F::Ring)
   iseven(n) || throw(ArgumentError("The dimension must be even"))
   G = MatrixGroup(n,F)
   G.descr = :Sp
   return G
end

function symplectic_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return symplectic_group(n, GF(b,a)[1])
end

function symplectic_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.Sp(_gap_filter(T), n, q))
   elseif T<:MatrixGroup   return symplectic_group(n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

function symplectic_group(::Type{T}, n::Int, F::Ring) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.Sp(_gap_filter(T), n, Int(order(F))))
   elseif T<:MatrixGroup   return symplectic_group(n,F)
   else throw(ArgumentError("Type not supported"))
   end
end

"""
    orthogonal_group(e::Int, n::Int, F::Ring)
    orthogonal_group(e::Int, n::Int, q::Int)
    GO = orthogonal_group

Return the orthogonal group of dimension `n` either over the field `F` or the field `GF(q)` of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted.
"""
function orthogonal_group(e::Int, n::Int, F::Ring)
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("GO+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("GO-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,F)
      G.descr = :GO
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

function orthogonal_group(e::Int, n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return orthogonal_group(e, n, GF(b,a)[1])
end

orthogonal_group(n::Int, F::Ring) = orthogonal_group(0,n,F)
orthogonal_group(n::Int, q::Int) = orthogonal_group(0,n,q)

function orthogonal_group(::Type{T}, e::Int, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.GO(_gap_filter(T), e, n, q))
   elseif T<:MatrixGroup   return orthogonal_group(e,n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

function orthogonal_group(::Type{T}, e::Int, n::Int, F::Ring) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.GO(_gap_filter(T), e, n, Int(order(F))))
   elseif T<:MatrixGroup   return orthogonal_group(e,n,F)
   else throw(ArgumentError("Type not supported"))
   end
end

orthogonal_group(::Type{T}, n::Int, F::Ring) where T<:GAPGroup = orthogonal_group(T,0,n,F)
orthogonal_group(::Type{T}, n::Int, q::Int) where T<:GAPGroup = orthogonal_group(T,0,n,q)

"""
    special_orthogonal_group(e::Int, n::Int, F::Ring)
    special_orthogonal_group(e::Int, n::Int, q::Int)
    SO = special_orthogonal_group

Return the special orthogonal group of dimension `n` either over the field `F` or the field `GF(q)` of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted.
"""
function special_orthogonal_group(e::Int, n::Int, F::Ring)
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("SO+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,F)
      G.descr = Symbol("SO-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,F)
      G.descr = :SO
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

function special_orthogonal_group(e::Int, n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return special_orthogonal_group(e, n, GF(b,a)[1])
end

special_orthogonal_group(n::Int, F::Ring) = special_orthogonal_group(0,n,F)
special_orthogonal_group(n::Int, q::Int) = special_orthogonal_group(0,n,q)

function special_orthogonal_group(::Type{T}, e::Int, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.SO(_gap_filter(T), e, n, q))
   elseif T<:MatrixGroup   return special_orthogonal_group(e,n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

function special_orthogonal_group(::Type{T}, e::Int, n::Int, F::Ring) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.SO(_gap_filter(T), e, n, Int(order(F))))
   elseif T<:MatrixGroup   return special_orthogonal_group(e,n,F)
   else throw(ArgumentError("Type not supported"))
   end
end

special_orthogonal_group(::Type{T}, n::Int, F::Ring) where T<:GAPGroup = special_orthogonal_group(T,0,n,F)
special_orthogonal_group(::Type{T}, n::Int, q::Int) where T<:GAPGroup = special_orthogonal_group(T,0,n,q)

"""
    omega_group(e::Int, n::Int, q::Int)

Return the Omega group of dimension `n` over the field `GF(q)` of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted.
"""
function omega_group(e::Int, n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   if e==1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,GF(b,a)[1])
      G.descr = Symbol("Omega+")
   elseif e==-1
      iseven(n) || throw(ArgumentError("The dimension must be even"))
      G = MatrixGroup(n,GF(b,a)[1])
      G.descr = Symbol("Omega-")
   elseif e==0
      isodd(n) || throw(ArgumentError("The dimension must be odd"))
      G = MatrixGroup(n,GF(b,a)[1])
      G.descr = :Omega
   else
      throw(ArgumentError("Invalid description of orthogonal group"))
   end
   return G
end

omega_group(n::Int, q::Int) = omega_group(0,n,q)

"""
    unitary_group(n::Int, q::Int)
    GU = unitary_group

Return the unitary group of dimension `n` over the field `GF(q^2)`.
"""
function unitary_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   G = MatrixGroup(n,GF(b,2*a)[1])
   G.descr = :GU
   return G
end

function unitary_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.GU(_gap_filter(T), n, q))
   elseif T<:MatrixGroup   return unitary_group(n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

"""
    special_unitary_group(n::Int, q::Int)
    SU = special_unitary_group

Return the special unitary group of dimension `n` over the field `GF(q^2)`.
"""
function special_unitary_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   G = MatrixGroup(n,GF(b,2*a)[1])
   G.descr = :SU
   return G
end

function special_unitary_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   if T<:PermGroup   return T(GAP.Globals.SU(_gap_filter(T), n, q))
   elseif T<:MatrixGroup   return special_unitary_group(n,q)
   else throw(ArgumentError("Type not supported"))
   end
end

const GL = general_linear_group
const SL = special_linear_group
const Sp = symplectic_group
const GO = orthogonal_group
const SO = special_orthogonal_group
const GU = unitary_group
const SU = special_unitary_group


"""
    matrix_group(V::T...) where T<:MatrixGroup
    matrix_group(V::AbstractVector{T}) where T<:MatrixGroup

Return the matrix group generated by elements in the vector `V`.
"""
function matrix_group(V::AbstractVector{T}) where T<:Union{MatElem,MatrixGroupElem}
   F = base_ring(V[1])
   n = nrows(V[1])
   G = GL(n,F)
   L = Vector{MatrixGroupElem}(undef, length(V))
   for i in 1:length(V)
      @assert base_ring(V[i])==F "The elements must have the same base ring"
      @assert nrows(V[1])==n "The elements must have the same dimension"
      if T<:MatElem
         L[i] = MatrixGroupElem(G,V[i])
      else
         L[i] = MatrixGroupElem(G,V[i].elm)
      end
   end
   H = MatrixGroup(n,F,L)
   for i in 1:length(V)   H.gens[i].parent = H   end
   return H
end

matrix_group(V::T...) where T<:Union{MatElem,MatrixGroupElem} = matrix_group(collect(V))




########################################################################
#
# Subgroups
#
########################################################################

function sub(G::MatrixGroup, elements::Vector{S}) where S <: GAPGroupElem
   @assert elem_type(G) == S
   elems_in_GAP = GAP.julia_to_gap(GapObj[x.X for x in elements])
   H = GAP.Globals.Group(elems_in_GAP)
   #H is the group. I need to return the inclusion map too
   K,f = _as_subgroup(G, H)
   L = Vector{MatrixGroupElem}(undef, length(elements))
   for i in 1:length(L)
      L[i] = MatrixGroupElem(K, elements[i].elm, elements[i].X)
   end
   K.gens = L
   return K,f
end



########################################################################
#
# Conjugation
#
########################################################################

function Base.show(io::IO, x::GroupConjClass{T,S}) where T <: MatrixGroup where S <: MatrixGroupElem
  show(io, x.repr)
  print(" ^ ")
  show(io, x.X)
end

function Base.show(io::IO, x::GroupConjClass{T,S}) where T <: MatrixGroup where S <: MatrixGroup
  show(io, x.repr)
  print(" ^ ")
  show(io, x.X)
end

function Base.:^(H::MatrixGroup, y::MatrixGroupElem)
   if isdefined(H,:gens)
      K = MatrixGroup(H.deg, H.ring)
      K.gens = [y^-1*x*y for x in H.gens]
      return K
   else
      return MatrixGroup(H.deg,H.ring,H.X^y.X)
   end
end

function conjugacy_classes_subgroups(G::MatrixGroup)
   L=GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,MatrixGroup(G.deg,G.ring,GAP.Globals.Representative(cc)),cc) for cc in L]
end

function conjugacy_classes_maximal_subgroups(G::MatrixGroup)
   L = GAP.gap_to_julia(Vector{GapObj},GAP.Globals.ConjugacyClassesMaximalSubgroups(G.X))
   return GroupConjClass{typeof(G), typeof(G)}[ _conjugacy_class(G,MatrixGroup(G.deg,G.ring,GAP.Globals.Representative(cc)),cc) for cc in L]
end

function Base.rand(C::GroupConjClass{S,T}) where S<:MatrixGroup where T<:MatrixGroup
   return MatrixGroup(C.X.deg,C.X.ring,GAP.Globals.Random(C.CC))
end
