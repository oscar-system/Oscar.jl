## Lazy

abstract type Lazy{T} end


### Const

struct Const{T} <: Lazy{T}
    c::T
end

Base.show(io::IO, c::Const) = print(io, c.c)


### Gen

struct Gen{T} <: Lazy{T}
    g::Symbol
end

Base.show(io::IO, g::Gen) = print(io, g.g)


### Plus

struct Plus{T} <: Lazy{T}
    xs::Vector{Lazy{T}}
end

Plus(xs::Lazy{T}...) where {T} = Plus(collect(Lazy{T}, xs))

function Base.show(io::IO, p::Plus)
    print(io, '(')
    join(io, p.xs, " + ")
    print(io, ')')
end


### Minus

struct Minus{T} <: Lazy{T}
    p::Lazy{T}
    q::Lazy{T}
end

Base.show(io::IO, p::Minus) = print(io, '(', p.p, " - ", p.q, ')')


### UniMinus

struct UniMinus{T} <: Lazy{T}
    p::Lazy{T}
end

Base.show(io::IO, p::UniMinus) = print(io, "(-", p.p, ')')


### Times

struct Times{T} <: Lazy{T}
    xs::Vector{Lazy{T}}
end

Times(xs::Lazy{T}...) where {T} = Times(collect(Lazy{T}, xs))

function Base.show(io::IO, p::Times)
    print(io, '(')
    join(io, p.xs)
    print(io, ')')
end


### Exp

struct Exp{T} <: Lazy{T}
    p::Lazy{T}
    e::Int
end

Base.show(io::IO, p::Exp) = print(io, p.p, '^', p.e)


### binary ops

#### +

+(x::Lazy{T}, y::Lazy{T}) where {T} = Plus(x, y)

function +(x::Plus{T}, y::Lazy{T}) where {T}
    p = Plus(copy(x.xs))
    push!(p.xs, y)
    p
end

function +(x::Lazy{T}, y::Plus{T}) where {T}
    p = Plus(copy(y.xs))
    pushfirst!(p.xs, x)
    p
end

function +(x::Plus{T}, y::Plus{T}) where {T}
    p = Plus(copy(x.xs))
    append!(p.xs, y.xs)
    p
end


#### -

-(p::Lazy{T}, q::Lazy{T}) where {T} = Minus(p, q)
-(p::Lazy{T}) where {T} = UniMinus(p)


#### *

*(x::Lazy{T}, y::Lazy{T}) where {T} = Times(x, y)

function *(x::Times{T}, y::Lazy{T}) where {T}
    p = Times(copy(x.xs))
    push!(p.xs, y)
    p
end

function *(x::Lazy{T}, y::Times{T}) where {T}
    p = Times(copy(y.xs))
    pushfirst!(p.xs, x)
    p
end

function *(x::Times{T}, y::Times{T}) where {T}
    p = Times(copy(x.xs))
    append!(p.xs, y.xs)
    p
end


#### ^

^(x::Lazy, e::Integer) = Exp(x, Int(e))


### adhoc binary ops

*(x, y::Lazy{T}) where {T} = Const(convert(T, x)) * y
*(x::Lazy{T}, y) where {T} = x * Const(convert(T, y))
