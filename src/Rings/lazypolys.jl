using .StraightLinePrograms: SLProgram
const SLP = StraightLinePrograms

## LazyPolyRing

struct LazyPolyRing{T<:RingElement,R<:Ring} <: MPolyRing{T}
    base_ring::R

    LazyPolyRing(r::Ring) = new{elem_type(r),typeof(r)}(r)
end

base_ring(F::LazyPolyRing) = F.base_ring

base_ring_type(::Type{LazyPolyRing{T, R}}) where {T,R} = R

## LazyPoly

struct LazyPoly{T<:RingElement,PR<:MPolyRing{T}} <: MPolyRingElem{T}
    parent::PR
    p::SLP.LazyRec
end

parent(p::LazyPoly) = p.parent

gen(R::LazyPolyRing, s::Symbol) = LazyPoly(R, SLP.Gen(s))

(R::LazyPolyRing)(s::Symbol) = gen(R, s)

(R::LazyPolyRing{T})(c::T) where {T} = LazyPoly(R, SLP.Const(c))

function check_parent(p::LazyPoly, q::LazyPoly)
    par = parent(p)
    @req par === parent(q) "incompatible parents"
    par
end

Base.show(io::IO, p::LazyPoly) = show(io, p.p)


### ops

+(p::LazyPoly, q::LazyPoly) = LazyPoly(check_parent(p, q), p.p + q.p)
-(p::LazyPoly, q::LazyPoly) = LazyPoly(check_parent(p, q), p.p - q.p)
-(p::LazyPoly) = LazyPoly(parent(p), -p.p)
*(p::LazyPoly, q::LazyPoly) = LazyPoly(check_parent(p, q), p.p * q.p)
^(p::LazyPoly, e::Base.Integer) = LazyPoly(parent(p), p.p^e)
