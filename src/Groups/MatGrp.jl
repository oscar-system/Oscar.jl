import AbstractAlgebra: MatElem, Ring
import Hecke: base_ring, det, fq_nmod, FqNmodFiniteField, tr, trace
import GAP: FFE

export
    GroupMatrix,
    GroupMatrixElem

mutable struct GroupMatrix{T<:Ring}
   X::GapObj
   gens::Vector{MatElem}                      # TODO: hard to put MatElem{T}; often, T does not coincide with typeof(ring)
   ring::T
   deg::Int
   descr::Symbol                       # e.g. GL
   ring_iso::Map{T,GapObj,SetMap}
   mat_iso::Map{MatElem{T},GapObj,SetMap}

   function GroupMatrix{T}(m::Int64, F::T) where T
      r = new{T}()
      r.deg = m
      r.ring = F
      return r
   end
end

mutable struct GroupMatrixElem{T<:Ring}
   parent::GroupMatrix{T}
   elm::MatElem
   X::GapObj

   function GroupMatrixElem{T}(G::GroupMatrix{T}, x::MatElem{S}) where {S,T}
      z = new{T}()
      z.parent = G
      z.elm = x
      return z
   end

   function GroupMatrixElem{T}(G::GroupMatrix{T}, x::MatElem{S}, x_gap::GapObj) where {S,T}
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

function Base.show(io::IO, x::GroupMatrix)
   if isdefined(x, :descr)
      print(io, string(x.descr), "(",x.deg,",",order(x.ring),")")
   else
      print(io, "Matrix group of degree ", x.deg, " over ", x.ring)
   end
end

function Base.show(io::IO, x::GroupMatrixElem)          #TODO: is this correct?
#   show(io)
   display(x.elm)
end

# add field :X
function get_gap_group(G::GroupMatrix{T}) where T
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
function get_gens(G::GroupMatrix{T}) where T
   if !isdefined(G,:gens)
      get_gap_group(G)
      F_gap = GAP.Globals.FieldOfMatrixGroup(G.X)
      L = GAP.Globals.GeneratorsOfGroup(G.X)
      G.gens = [MatGapToOscar(L[i],F_gap,G.ring) for i in 1:length(L)]
   end
end

# add field :X
function get_gap_matrix(x::GroupMatrixElem{T}) where T
   if !isdefined(x,:X)
      x.X = MatOscarToGap(x.elm, x.parent.ring)
   end
end

Base.IteratorSize(::Type{<:GroupMatrix}) = Base.SizeUnknown()

function Base.iterate(G::GroupMatrix{T}) where T
  get_gap_group(G)
  L=GAP.Globals.Iterator(G.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return GroupMatrixElem{T}(G, MatGapToOscar(i,GAP.Globals.DefaultFieldOfMatrix(i),G.ring),i), L
end

function Base.iterate(G::GroupMatrix{T}, state) where T
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return GroupMatrixElem{T}(G, MatGapToOscar(i,GAP.Globals.DefaultFieldOfMatrix(i),G.ring),i), state
end

# need this function just for the iterator
Base.length(x::GroupMatrix)::Int = order(x)

########################################################################
#
# Membership
#
########################################################################

# this saves the value of x.X
# x_gap = x.X if this is already known, x_gap = 0 otherwise
function lies_in(x::MatElem, G::GroupMatrix, x_gap)
   if x==identity_matrix(G.ring,G.deg) return true, Nothing end
   if isdefined(G,:gens)
      if x in G.gens return true, Nothing end
   else
      get_gap_group(G)
      if x_gap==0 x_gap = MatOscarToGap(x,base_ring(x)) end
     # x_gap !=0 || x_gap = MatOscarToGap(x,base_ring(x))
      return GAP.Globals.in(x_gap,G.X), x_gap
   end
end

Base.in(x::MatElem, G::GroupMatrix) = lies_in(x,G,0)[1]

function Base.in(x::GroupMatrixElem, G::GroupMatrix)
   if isdefined(x,:X) return lies_in(x.elm,G,x.X)[1]
   else return lies_in(x.elm,G,0)[1]
   end
end

# embedding an element of type MatElem into a group G
function (G::GroupMatrix{T})(x::MatElem) where T
   vero, x_gap = lies_in(x,G,0)
   vero || throw(ArgumentError("Element not in the group"))
   if x_gap==Nothing return GroupMatrixElem{T}(G,x)
   else return GroupMatrixElem{T}(G,x,x_gap)
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
# Existing functions on elements
#
########################################################################

det(x::GroupMatrixElem{T}) where T = det(x.elm)

trace(x::GroupMatrixElem{T}) where T = trace(x.elm)
tr(x::GroupMatrixElem{T}) where T = tr(x.elm)

function order(x::GroupMatrixElem{T}) where T
   get_gap_matrix(x)
   return GAP.Globals.Order(x.X)
end

########################################################################
#
# Methods on groups
#
########################################################################

base_ring(G::GroupMatrix{T}) where T = G.ring

function Base.rand(G::GroupMatrix{T}) where T
   get_gap_group(G)
   x_gap = GAP.Globals.Random(G.X)
   x_oscar = MatGapToOscar(x_gap,GAP.Globals.FieldOfMatrixGroup(G.X),G.ring)
   x = GroupMatrixElem{T}(G,x_oscar,x_gap)
   return x
end

function gens(G::GroupMatrix{T}) where T
   get_gens(G)
   V = [GroupMatrixElem{T}(G,g) for g in G.gens]
   return V
end

gens(G::GroupMatrix{T}, i::Int) where T = gens(G)[i]

Base.getindex(G::GroupMatrix{T}, i::Int) where T = gens(G, i)

function ngens(G::GroupMatrix{T}) where T
   get_gens(G)
   return length(G.gens)
end

function order(G::GroupMatrix)
   get_gap_group(G)
   return GAP.Globals.Order(G.X)
end

########################################################################
#
# Constructors
#
########################################################################

function general_linear_group(n::Int, F::Ring)
   G = GroupMatrix{typeof(F)}(n,F)
   G.descr = :GL
   return G
end

function special_linear_group(n::Int, F::Ring)
   G = GroupMatrix{typeof(F)}(n,F)
   G.descr = :SL
   return G
end


