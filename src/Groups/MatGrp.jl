import AbstractAlgebra: MatElem, Ring
import Hecke: base_ring, FqNmodFiniteField

export
    GroupMatrix,
    GroupMatrixElem

mutable struct GroupMatrix{T<:Ring}
   X::GapObj
   gens::Vector{MatElem{T}}
   ring::T
   deg::Int
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
end

Base.show(io::IO, x::GroupMatrix) = print(io, "Matrix group of degree ", x.deg, " over ", x.ring)
Base.show(io::IO, x::GroupMatrixElem) = print(io, x.elm)



########################################################################
#
# conversions
#
########################################################################

function FieldGapToJulia(F::GapObj)
   p = GAP.Globals.Characteristic(F)
   q = GAP.Globals.Size(F)
   d = Base.Integer(log(p,q))

   return GF(p,d)
end

function FieldJuliaToGap(F::FqNmodFiniteField)
   p = Int64(characteristic(F))
   d = Int64(degree(F))

   return GAP.Globals.GF(p,d)
end


########################################################################
#
# Methods
#
########################################################################

base_ring(G::GroupMatrix{T}) where T = G.ring
