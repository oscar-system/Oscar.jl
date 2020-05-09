## SLPolyRing (SL = straight-line)

struct SLPolyRing{T<:RingElement,R<:Ring} <: MPolyRing{T}
    base_ring::R
    S::Vector{Symbol}

    SLPolyRing(r::Ring, s::Vector{Symbol}) = new{elem_type(r),typeof(r)}(r, s)
end

base_ring(S::SLPolyRing) = S.base_ring

symbols(S::SLPolyRing) = S.S

(PR::SLPolyRing{T})() where {T} = SLPoly(PR, T[], UInt64[])


## SLPoly

struct SLPoly{T<:RingElement,SLPR<:SLPolyRing{T}} <: MPolyElem{T}
    parent::SLPR
    cs::Vector{T}          # constants
    lines::Vector{UInt64}  # instructions
    f::Ref{Function}       # compiled evalutation

    SLPoly(parent, xs, c) =
        new{elem_type(base_ring(parent)),typeof(parent)}(
            parent, xs, c, Ref{Function}())
end

parent(p::SLPoly) = p.parent
