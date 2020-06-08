## Lazy

abstract type Lazy end

gens(l::Lazy) = sort!(unique!(pushgens!(Symbol[], l)::Vector{Symbol}))

Base.:(==)(k::Lazy, l::Lazy) = false

Lazy(x::Lazy) = x
Lazy(x) = Const(x)

evaluate(l::Lazy, xs) = evaluate(gens(l), l, xs)


### Const

struct Const{T} <: Lazy
    c::T
end

Base.show(io::IO, c::Const) = print(io, c.c)

pushgens!(gs, c::Const) = gs
constantstype(l::Const{T}) where {T} = T

Base.:(==)(k::Const, l::Const) = k.c == l.c

evaluate(gs, c::Const, xs) = c.c


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

evaluate(gs, g::Gen, xs) = xs[findfirst(==(g.g), gs)]


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

evaluate(gs, p::Plus, xs) = mapreduce(q -> evaluate(gs, q, xs), +, p.xs)


### Minus

struct Minus <: Lazy
    p::Lazy
    q::Lazy
end

Base.show(io::IO, p::Minus) = print(io, '(', p.p, " - ", p.q, ')')

pushgens!(gs, l::Minus) = pushgens!(pushgens!(gs, l.p), l.q)

constantstype(m::Minus) = typejoin(constantstype(m.p), constantstype(m.q))

Base.:(==)(k::Minus, l::Minus) = k.p == l.p && k.q == l.q

evaluate(gs, p::Minus, xs) = evaluate(gs, p.p, xs) - evaluate(gs, p.q, xs)


### UniMinus

struct UniMinus <: Lazy
    p::Lazy
end

Base.show(io::IO, p::UniMinus) = print(io, "(-", p.p, ')')

pushgens!(gs, l::UniMinus) = pushgens!(gs, l.p)

constantstype(p::UniMinus) = constantstype(p.p)

Base.:(==)(k::UniMinus, l::UniMinus) = k.p == l.p

evaluate(gs, p::UniMinus, xs) = -evaluate(gs, p.p, xs)


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

evaluate(gs, p::Times, xs) = mapreduce(q -> evaluate(gs, q, xs), *, p.xs)


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

evaluate(gs, p::Exp, xs) = evaluate(gs, p.p, xs)^p.e


### Decision

struct Decision <: Lazy
    ps::Vector{Tuple{Lazy,Int}}
end

test(p::Lazy, n) = Decision(Tuple{Lazy,Int}[(p, n)])

Base.show(io::IO, d::Decision) =
    join(io, ("test($p, $n)" for (p, n) in d.ps), " & ")

constantstype(l::Decision) =
    mapreduce(constantstype∘first, typejoin, l.ps, init=Union{})

pushgens!(gs, l::Decision) =
    foldl((gs, d) -> pushgens!(gs, first(d)), l.ps, init=gs)

function Base.:(&)(p::Decision, q::Decision)
    d = Decision(copy(p.ps))
    append!(d.ps, q.ps)
    d
end

Base.:(==)(p::Decision, q::Decision) = p.ps == q.ps

function evaluate(gs, p::Decision, xs)
    (d, i), rest = Iterators.peel(p.ps) # p.ps should never be empty
    res = test(evaluate(gs, d, xs), i)
    res === false && return false
    for (d, i) in rest
        r = test(evaluate(gs, d, xs), i)
        r === false && return false
        res &= r
    end
    res
end


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
