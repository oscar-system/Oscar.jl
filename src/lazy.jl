## Lazy

abstract type Lazy end


### Const

struct Const{T} <: Lazy
    c::T
end

Base.show(io::IO, c::Const) = print(io, c.c)


### Gen

struct Gen <: Lazy
    g::Symbol
end

Base.show(io::IO, g::Gen) = print(io, g.g)


### Plus

struct Plus <: Lazy
    xs::Vector{Lazy}
end

Plus(xs::Lazy...) = Plus(collect(Lazy, xs))

function Base.show(io::IO, p::Plus)
    print(io, '(')
    join(io, p.xs, " + ")
    print(io, ')')
end


### Minus

struct Minus <: Lazy
    p::Lazy
    q::Lazy
end

Base.show(io::IO, p::Minus) = print(io, '(', p.p, " - ", p.q, ')')


### UniMinus

struct UniMinus <: Lazy
    p::Lazy
end

Base.show(io::IO, p::UniMinus) = print(io, "(-", p.p, ')')


### Times

struct Times <: Lazy
    xs::Vector{Lazy}
end

Times(xs::Lazy...) = Times(collect(Lazy, xs))

function Base.show(io::IO, p::Times)
    print(io, '(')
    join(io, p.xs)
    print(io, ')')
end


### Exp

struct Exp <: Lazy
    p::Lazy
    e::Int
end

Base.show(io::IO, p::Exp) = print(io, p.p, '^', p.e)


### binary ops

#### +

+(x::Lazy, y::Lazy) = Plus(x, y)

function +(x::Plus, y::Lazy)
    p = Plus(copy(x.xs))
    push!(p.xs, y)
    p
end

function +(x::Lazy, y::Plus)
    p = Plus(copy(y.xs))
    pushfirst!(p.xs, x)
    p
end

function +(x::Plus, y::Plus)
    p = Plus(copy(x.xs))
    append!(p.xs, y.xs)
    p
end


#### -

-(p::Lazy, q::Lazy) = Minus(p, q)
-(p::Lazy) = UniMinus(p)


#### *

*(x::Lazy, y::Lazy) = Times(x, y)

function *(x::Times, y::Lazy)
    p = Times(copy(x.xs))
    push!(p.xs, y)
    p
end

function *(x::Lazy, y::Times)
    p = Times(copy(y.xs))
    pushfirst!(p.xs, x)
    p
end

function *(x::Times, y::Times)
    p = Times(copy(x.xs))
    append!(p.xs, y.xs)
    p
end


#### ^

^(x::Lazy, e::Integer) = Exp(x, Int(e))


### adhoc binary ops

*(x, y::Lazy) = Const(x) * y
*(x::Lazy, y) = x * Const(y)
