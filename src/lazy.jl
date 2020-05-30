## Lazy

abstract type Lazy end

gens(l::Lazy) = sort!(unique!(pushgens!(Symbol[], l)::Vector{Symbol}))

Base.:(==)(k::Lazy, l::Lazy) = false

### Const

struct Const{T} <: Lazy
    c::T
end

Base.show(io::IO, c::Const) = print(io, c.c)

pushgens!(gs, c::Const) = gs
constantstype(l::Const{T}) where {T} = T

Base.:(==)(k::Const, l::Const) = k.c == l.c


### Gen

struct Gen <: Lazy
    g::Symbol
end

Base.show(io::IO, g::Gen) = print(io, g.g)

pushgens!(gs, g::Gen) =
    if findfirst(==(g.g), gs) === nothing
        push!(gs, g.g)
    else
        gs
    end

constantstype(l::Gen) = Union{}

Base.:(==)(k::Gen, l::Gen) = k.g == l.g


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

pushgens!(gs, p::Plus) = foldl(pushgens!, p.xs, init=gs)

constantstype(p::Plus) =
    mapreduce(constantstype, typejoin, p.xs, init=Union{})

Base.:(==)(k::Plus, l::Plus) = k.xs == l.xs


### Minus

struct Minus <: Lazy
    p::Lazy
    q::Lazy
end

Base.show(io::IO, p::Minus) = print(io, '(', p.p, " - ", p.q, ')')

pushgens!(gs, l::Minus) = pushgens!(pushgens!(gs, l.p), l.q)

constantstype(m::Minus) = typejoin(constantstype(m.p), constantstype(m.q))

Base.:(==)(k::Minus, l::Minus) = k.p == l.p && k.q == l.q


### UniMinus

struct UniMinus <: Lazy
    p::Lazy
end

Base.show(io::IO, p::UniMinus) = print(io, "(-", p.p, ')')

pushgens!(gs, l::UniMinus) = pushgens!(gs, l.p)

constantstype(p::UniMinus) = constantstype(p.p)

Base.:(==)(k::UniMinus, l::UniMinus) = k.p == l.p


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

pushgens!(gs, p::Times) = foldl(pushgens!, p.xs, init=gs)

constantstype(p::Times) =
    mapreduce(constantstype, typejoin, p.xs, init=Union{})

Base.:(==)(k::Times, l::Times) = k.xs == l.xs


### Exp

struct Exp <: Lazy
    p::Lazy
    e::Int
end

Base.show(io::IO, p::Exp) = print(io, p.p, '^', p.e)

pushgens!(gs, l::Exp) = pushgens!(gs, l.p)

constantstype(p::Exp) = constantstype(p.p)

Base.:(==)(k::Exp, l::Exp) = k.p == l.p && k.e == l.e

Base.literal_pow(::typeof(^), p::Lazy, ::Val{e}) where {e} = Exp(p, e)

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

+(x, y::Lazy) = Const(x) + y
+(x::Lazy, y) = x + Const(y)

-(x, y::Lazy) = Const(x) - y
-(x::Lazy, y) = x - Const(y)
