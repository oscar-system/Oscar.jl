import AbstractAlgebra: MatElem, MatSpace, Ring
import Hecke: base_ring, det, fq_nmod, FqNmodFiniteField, tr, trace
import GAP: FFE

export
    MatrixGroup,
    MatrixGroupElem

# isomorphism from the Oscar field to the GAP field
mutable struct GenRingIso{S<:Ring,T<:GapObj} <: Map{S,T,SetMap,GenRingIso}
   domain::S
   codomain::T
   f
end

# sends the Oscar matrix into the GAP matrix
mutable struct GenMatIso{S<:MatSpace, T<:GapObj} <: Map{S,T,SetMap,GenMatIso{S,T}}
   domain::S
   codomain::T
   f
end

mutable struct MatrixGroup{T<:Ring}
   X::GapObj
   gens::Vector                     # TODO: hard to put MatElem{T}; often, T does not coincide with typeof(ring)
   ring::T
   deg::Int
   descr::Symbol                       # e.g. GL
   ring_iso::Map{T,GapObj,SetMap}
   mat_iso::Map{MatElem{T},GapObj,SetMap}

   function MatrixGroup{T}(m::Int64, F::T) where T
      r = new{T}()
      r.deg = m
      r.ring = F
      return r
   end
end

mutable struct MatrixGroupElem{T<:Ring}
   parent::MatrixGroup{T}
   elm::MatElem
   X::GapObj

   function MatrixGroupElem{T}(G::MatrixGroup{T}, x::MatElem{S}) where {S,T}
      z = new{T}()
      z.parent = G
      z.elm = x
      return z
   end

   function MatrixGroupElem{T}(G::MatrixGroup{T}, x::MatElem{S}, x_gap::GapObj) where {S,T}
      z = new{T}()
      z.parent = G
      z.elm = x
      z.X = x_gap
      return z
   end
end



########################################################################
#
# Basic
#
########################################################################

function Base.show(io::IO, x::MatrixGroup)
   if isdefined(x, :descr)
      print(io, string(x.descr), "(",x.deg,",",order(x.ring),")")
   else
      print(io, "Matrix group of degree ", x.deg, " over ", x.ring)
   end
end

function Base.show(io::IO, x::MatrixGroupElem)          #TODO: is this correct?
#   show(io)
   display(x.elm)
end

function Base.getproperty(G::MatrixGroup{T}, sym::Symbol) where T
   if sym === :X
      if isdefined(G,:X) return getfield(G,:X) end
      if isdefined(G,:descr)
         if G.descr==:GL setfield!(G,:X,GAP.Globals.GL(G.deg,Int(order(G.ring))))
         elseif G.descr==:SL setfield!(G,:X,GAP.Globals.SL(G.deg,Int(order(G.ring))))
         end
      elseif isdefined(G,:gens)
         V = GAP.julia_to_gap([MatOscarToGap(g,g.ring) for g in G.gens])
         setfield!(G,:X,GAP.Globals.Group(V))
      end
   elseif sym === :gens
      if isdefined(G, :gens) return getfield(G, :gens) end
      if !isdefined(G, :X)
         if G.descr==:GL setfield!(G,:X,GAP.Globals.GL(G.deg,Int(order(G.ring))))
         elseif G.descr==:SL setfield!(G,:X,GAP.Globals.SL(G.deg,Int(order(G.ring))))
         end
      end
      F_gap = GAP.Globals.FieldOfMatrixGroup(getfield(G,:X))
      L = GAP.Globals.GeneratorsOfGroup(getfield(G,:X))
      setfield!(G,:gens,[MatGapToOscar(L[i],F_gap,G.ring) for i in 1:length(L)])
   else
      return getfield(G, sym)
   end
end

#=
# add field :X
function get_gap_group(G::MatrixGroup{T}) where T
   if !isdefined(G,:X)
      if isdefined(G,:descr)
         if G.descr==:GL G.X = GAP.Globals.GL(G.deg,Int(order(G.ring)))
         elseif G.descr==:SL G.X = GAP.Globals.SL(G.deg,Int(order(G.ring)))
         end
      elseif isdefined(G,:gens)
         V = GAP.julia_to_gap([MatOscarToGap(g,g.ring) for g in G.gens])
         G.X = GAP.Globals.Group(V)
      end
   end
end

# add field :gens
function get_gens(G::MatrixGroup{T}) where T
   if !isdefined(G,:gens)
      get_gap_group(G)
      F_gap = GAP.Globals.FieldOfMatrixGroup(G.X)
      L = GAP.Globals.GeneratorsOfGroup(G.X)
      G.gens = [MatGapToOscar(L[i],F_gap,G.ring) for i in 1:length(L)]
   end
end
=#

# add field :X
function Base.getproperty(x::MatrixGroupElem{T}, sym::Symbol) where T
   if sym === :X && !isdefined(x,:X)
      setfield!(x, :X, MatOscarToGap(x.elm, x.parent.ring))
   end
   return getfield(x,sym)
end

Base.IteratorSize(::Type{<:MatrixGroup}) = Base.SizeUnknown()

function Base.iterate(G::MatrixGroup{T}) where T
  L=GAP.Globals.Iterator(G.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return MatrixGroupElem{T}(G, MatGapToOscar(i,GAP.Globals.DefaultFieldOfMatrix(i),G.ring),i), L
end

function Base.iterate(G::MatrixGroup{T}, state) where T
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return MatrixGroupElem{T}(G, MatGapToOscar(i,GAP.Globals.DefaultFieldOfMatrix(i),G.ring),i), state
end

# need this function just for the iterator
Base.length(x::MatrixGroup)::Int = order(x)

########################################################################
#
# Membership
#
########################################################################

# this saves the value of x.X
# x_gap = x.X if this is already known, x_gap = 0 otherwise
function lies_in(x::MatElem, G::MatrixGroup, x_gap)
   if isone(x) return true, Nothing end
   if isdefined(G,:gens) && x in G.gens
      return true, Nothing
   else
      if x_gap==0 x_gap = MatOscarToGap(x,base_ring(x)) end
     # x_gap !=0 || x_gap = MatOscarToGap(x,base_ring(x))
      return GAP.Globals.in(x_gap,G.X), x_gap
   end
end

Base.in(x::MatElem, G::MatrixGroup) = lies_in(x,G,0)[1]

function Base.in(x::MatrixGroupElem, G::MatrixGroup)
   if isdefined(x,:X) return lies_in(x.elm,G,x.X)[1]
   else return lies_in(x.elm,G,0)[1]
   end
end

# embedding an element of type MatElem into a group G
function (G::MatrixGroup{T})(x::MatElem) where T
   vero, x_gap = lies_in(x,G,0)
   vero || throw(ArgumentError("Element not in the group"))
   if x_gap==Nothing return MatrixGroupElem{T}(G,x)
   else return MatrixGroupElem{T}(G,x,x_gap)
   end
end

########################################################################
#
# conversions
#
########################################################################
#TODO: support other field types

function FieldGapToOscar(F::GapObj)
   p = GAP.Globals.Characteristic(F)
   q = GAP.Globals.Size(F)
   d = Base.Integer(log(p,q))

   return GF(p,d)
end

function FieldOscarToGap(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))

   return GAP.Globals.GF(p,d)
end

function FieldElemGapToOscar(x::FFE, Fgap::GapObj, F::FqNmodFiniteField)
   q = GAP.Globals.Size(Fgap)
   z = gen(F)
   if GAP.Globals.IsZero(x)
      return 0*z
   else
      d = Integer(GAP.Globals.LogFFE(x, GAP.Globals.Z(q)))
      return z^d
   end
end

function FieldElemGapToOscar(x::FFE, F::GapObj)
   q = GAP.Globals.Size(F)
   K,z = FieldGapToOscar(F)
   if GAP.Globals.IsZero(x)
      return 0*z
   else
      d = Integer(GAP.Globals.LogFFE(x, GAP.Globals.Z(q)))
      return z^d
   end
end

function FieldElemOscarToGap(x::fq_nmod, F::FqNmodFiniteField)
   q = Int64(order(F))
   d = degree(F)
   z = gen(F)
   v = [Int64(coeff(x,i)) for i in 0:d]
   y = sum([GAP.Globals.Z(q)^i*v[i+1] for i in 0:d])
   
   return y
end

function MatGapToOscar(x::GapObj, Fgap::GapObj, F::FqNmodFiniteField)
   n = GAP.Globals.Size(x)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [FieldElemGapToOscar(Arr[i][j],Fgap,F) for i in 1:n for j in 1:n]

   return matrix(F,n,n,L)
end

function MatGapToOscar(x::GapObj, F::GapObj)
   n = GAP.Globals.Size(x)
   Arr = [GAP.gap_to_julia(x[i]) for i in 1:n]
   L = [FieldElemGapToOscar(Arr[i][j],F) for i in 1:n for j in 1:n]

   return matrix(FieldGapToOscar(F)[1],n,n,L)
end

function MatOscarToGap(x::MatElem, F::FqNmodFiniteField)
   r = nrows(x)
   c = ncols(x)
   S = Vector{GapObj}(undef,r)
   for i in 1:r
      S[i] = GAP.julia_to_gap([FieldElemOscarToGap(x[i,j],F) for j in 1:c])
   end

   return GAP.julia_to_gap(S)
end

########################################################################
#
# Methods on elements
#
########################################################################

# TODO: we are not currently keeping track of the parent
function ==(x::MatrixGroupElem{T},y::MatrixGroupElem{T}) where T
   return x.elm==y.elm
end

# TODO: we wish to multiply / conjugate also matrices with different parents ?

function Base.:*(x::MatrixGroupElem{T},y::MatrixGroupElem{T}) where T
   if x.parent==y.parent
      return MatrixGroupElem{T}(x.parent, x.elm*y.elm)
   else
      throw(ArgumentError("Matrices not in the same group"))
   end
end

Base.:^(x::MatrixGroupElem{T}, n::Int) where T = MatrixGroupElem{T}(x.parent, x.elm^n)

Base.isone(x::MatrixGroupElem{T}) where T = isone(x.elm)

Base.inv(x::MatrixGroupElem{T}) where T = MatrixGroupElem{T}(x.parent, inv(x.elm))

function Base.conj(x::MatrixGroupElem{T}, y::MatrixGroupElem{T}) where T
   if x.parent==y.parent
      return MatrixGroupElem{T}(x.parent,inv(y.elm)*x.elm*y.elm)
   else
      throw(ArgumentError("Matrices not in the same group"))
   end
end

Base.:^(x::MatrixGroupElem{T}, y::MatrixGroupElem{T}) where T = conj(x,y)

comm(x::MatrixGroupElem{T}, y::MatrixGroupElem{T}) where T = inv(x)*conj(x,y)

det(x::MatrixGroupElem{T}) where T = det(x.elm)
base_ring(x::MatrixGroupElem) = x.parent.ring

parent(x::MatrixGroupElem{T}) where T = x.parent

Base.getindex(x::MatrixGroupElem{T}, i::Int, j::Int) where T = x.elm[i,j]

trace(x::MatrixGroupElem{T}) where T = trace(x.elm)
tr(x::MatrixGroupElem{T}) where T = tr(x.elm)

order(x::MatrixGroupElem{T}) where T = GAP.Globals.Order(x.X)

########################################################################
#
# Methods on groups
#
########################################################################

base_ring(G::MatrixGroup{T}) where T = G.ring

Base.one(G::MatrixGroup{T}) where T = MatrixGroupElem{T}(G, identity_matrix(G.ring, G.deg))

function Base.rand(G::MatrixGroup{T}) where T
   x_gap = GAP.Globals.Random(G.X)
   x_oscar = MatGapToOscar(x_gap,GAP.Globals.FieldOfMatrixGroup(G.X),G.ring)
   x = MatrixGroupElem{T}(G,x_oscar,x_gap)
   return x
end

function gens(G::MatrixGroup{T}) where T
   V = [MatrixGroupElem{T}(G,g) for g in G.gens]
   return V
end

gens(G::MatrixGroup{T}, i::Int) where T = gens(G)[i]

Base.getindex(G::MatrixGroup{T}, i::Int) where T = gens(G, i)

function ngens(G::MatrixGroup{T}) where T
   return length(G.gens)
end

function order(G::MatrixGroup)
   return GAP.Globals.Order(G.X)
end

########################################################################
#
# Constructors
#
########################################################################

function general_linear_group(n::Int, F::Ring)
   G = MatrixGroup{typeof(F)}(n,F)
   G.descr = :GL
   return G
end

function general_linear_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return general_linear_group(n, GF(b,a)[1])
end

function special_linear_group(n::Int, F::Ring)
   G = MatrixGroup{typeof(F)}(n,F)
   G.descr = :SL
   return G
end

function special_linear_group(n::Int, q::Int)
   (a,b) = ispower(q)
   if !isprime(b) throw(ArgumentError("The field size must be a prime power")) end
   return special_linear_group(n, GF(b,a)[1])
end

const GL = general_linear_group
const SL = special_linear_group


