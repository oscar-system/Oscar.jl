## LazyPolyRing

struct LazyPolyRing{T<:RingElement,R<:Ring} <: MPolyRing{T}
   base_ring::R

   LazyPolyRing(r::Ring) = new{elem_type(r),typeof(r)}(r)
end

base_ring(F::LazyPolyRing) = F.base_ring


## RecPoly

abstract type RecPoly{T<:RingElement} end


### Const

struct Const{T} <: RecPoly{T}
   c::T
end

Base.show(io::IO, c::Const) = print(io, c.c)


### Gen

struct Gen{T} <: RecPoly{T}
   g::Symbol
end

Base.show(io::IO, g::Gen) = print(io, g.g)


### PlusPoly

struct PlusPoly{T} <: RecPoly{T}
   xs::Vector{RecPoly{T}}
end

PlusPoly(xs::RecPoly{T}...) where {T} = PlusPoly(collect(RecPoly{T}, xs))

function Base.show(io::IO, p::PlusPoly)
   print(io, '(')
   join(io, p.xs, " + ")
   print(io, ')')
end


### MinusPoly

struct MinusPoly{T} <: RecPoly{T}
   p::RecPoly{T}
   q::RecPoly{T}
end

Base.show(io::IO, p::MinusPoly) = print(io, '(', p.p, " - ", p.q, ')')


### UniMinus

struct UniMinusPoly{T} <: RecPoly{T}
   p::RecPoly{T}
end

Base.show(io::IO, p::UniMinusPoly) = print(io, "(-", p.p, ')')


### TimesPoly

struct TimesPoly{T} <: RecPoly{T}
   xs::Vector{RecPoly{T}}
end

TimesPoly(xs::RecPoly{T}...) where {T} = TimesPoly(collect(RecPoly{T}, xs))

function Base.show(io::IO, p::TimesPoly)
   print(io, '(')
   join(io, p.xs)
   print(io, ')')
end


### ExpPoly

struct ExpPoly{T} <: RecPoly{T}
   p::RecPoly{T}
   e::Int
end

Base.show(io::IO, p::ExpPoly) = print(io, p.p, '^', p.e)


### binary ops

#### +

+(x::RecPoly{T}, y::RecPoly{T}) where {T} = PlusPoly(x, y)

function +(x::PlusPoly{T}, y::RecPoly{T}) where {T}
   p = PlusPoly(copy(x.xs))
   push!(p.xs, y)
   p
end

function +(x::RecPoly{T}, y::PlusPoly{T}) where {T}
   p = PlusPoly(copy(y.xs))
   pushfirst!(p.xs, x)
   p
end

function +(x::PlusPoly{T}, y::PlusPoly{T}) where {T}
   p = PlusPoly(copy(x.xs))
   append!(p.xs, y.xs)
   p
end


#### -

-(p::RecPoly{T}, q::RecPoly{T}) where {T} = MinusPoly(p, q)
-(p::RecPoly{T}) where {T} = UniMinusPoly(p)


#### *

*(x::RecPoly{T}, y::RecPoly{T}) where {T} = TimesPoly(x, y)

function *(x::TimesPoly{T}, y::RecPoly{T}) where {T}
   p = TimesPoly(copy(x.xs))
   push!(p.xs, y)
   p
end

function *(x::RecPoly{T}, y::TimesPoly{T}) where {T}
   p = TimesPoly(copy(y.xs))
   pushfirst!(p.xs, x)
   p
end

function *(x::TimesPoly{T}, y::TimesPoly{T}) where {T}
   p = TimesPoly(copy(x.xs))
   append!(p.xs, y.xs)
   p
end


#### ^

^(x::RecPoly, e::Integer) = ExpPoly(x, Int(e))


### adhoc binary ops

*(x::RingElement, y::RecPoly{T}) where {T} = Const(convert(T, x)) * y
*(x::RecPoly{T}, y::RingElement) where {T} = x * Const(convert(T, y))


## LazyPoly

struct LazyPoly{T<:RingElement,PR<:MPolyRing{T}} <: MPolyElem{T}
   parent::PR
   p::RecPoly{T}
end

parent(p::LazyPoly) = p.parent

gen(R::LazyPolyRing, s::Symbol) = LazyPoly(R, Gen{elem_type(R.base_ring)}(s))

(R::LazyPolyRing)(s::Symbol) = gen(R, s)

(R::LazyPolyRing{T})(c::T) where {T} = LazyPoly(R, Const(c))

Base.show(io::IO, p::LazyPoly) = show(io, p.p)
