## LazyPolyRing

struct LazyPolyRing{T<:RingElement,R<:Ring} <: MPolyRing{T}
   base_ring::R

   LazyPolyRing(r::Ring) = new{elem_type(r),typeof(r)}(r)
end

base_ring(F::LazyPolyRing) = F.base_ring
