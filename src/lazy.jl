## Free

abstract type Lazy end

struct Free <: AbstractSLProgram
    x::Lazy # must not contain Gen
    gens::Vector{Symbol}
end

Free(x::Lazy) = Free(x, collect(Symbol, slpsyms(maxinput(x))))

Free(x) = Free(Const(x), Symbol[])

function freegens(n::Integer, syms=slpsyms(n))
    if !isa(syms, Vector{Symbol})
        syms = collect(Symbol, syms)
    end
    Free[Free(Input(i), syms) for i=1:n]
end

freegens(syms::AbstractVector{Symbol}) = freegens(length(syms), syms)

gens(::Type{Free}, n::Integer) = freegens(n)

ngens(f::Free) = length(f.gens)

gens(f::Free) = f.gens

Base.show(io::IO, f::Free) =
    show(IOContext(io, :slp_free_gens => gens(f)), f.x)

evaluate(f::Free, xs) = evaluate(gens(f), f.x, xs)

# check compatibility and return the biggest of the two gens arrays
function gens(x::Free, xs::Free...)
    gx = gens(x)
    for y in xs
        gy = gens(y)
        if length(gy) > length(gx)
            gx, gy = gy, gx
        end
        gy == view(gx, eachindex(gy)) || throw(ArgumentError(
            "incompatible symbols"))
    end
    gx
end

# TODO: must they also have the same gens?
Base.:(==)(x::Free, y::Free) = x.x == y.x

compile(::Type{SLProgram}, f::Free) =
    compile(SLProgram{constantstype(f.x)}, f)

function compile(::Type{SLProgram{T}}, f::Free) where T
    p = SLProgram{T}()
    i = pushlazy!(p, f.x, gens(f))
    pushfinalize!(p, i)
end

compile(f::Free) = compile(SLProgram, f)


### unary/binary ops

+(x::Free, y::Free) = Free(x.x + y.x, gens(x, y))
-(x::Free, y::Free) = Free(x.x - y.x, gens(x, y))
*(x::Free, y::Free) = Free(x.x * y.x, gens(x, y))

Base.:(&)(x::Free, y::Free) = Free(x.x & y.x, gens(x, y))

-(x::Free) = Free(-x.x, gens(x))
^(x::Free, e::Integer) = Free(x.x^e, gens(x))

Base.literal_pow(::typeof(^), p::Free, ::Val{e}) where {e} = p^e

test(x::Free, n) = Free(test(x.x, n), gens(x))

# TODO: don't splat
list(::Type{Free}, xs) = Free(List(Lazy[x.x for x in xs]), gens(xs...))

function compose(p::Free, q::Free; flatten=true)
    if flatten
        q.x isa List || throw(ArgumentError(
            "first argument must return a list"))
        gs = gens(q)
        evaluate(p, [Free(qi, gs) for qi in q.x.xs])::Free
    else
        Free(Compose(p.x, q.x), gens(q))
    end
end

Base.getindex(p::Free, is::Free...) =
    Free(Getindex(p.x, Lazy[i.x for i in is]), gens(p, is...))


#### adhoc

+(x::Free, y) = Free(x.x + y, gens(x))
+(x, y::Free) = Free(x + y.x, gens(y))

-(x::Free, y) = Free(x.x - y, gens(x))
-(x, y::Free) = Free(x - y.x, gens(y))

*(x::Free, y) = Free(x.x * y, gens(x))
*(x, y::Free) = Free(x * y.x, gens(y))

Base.getindex(x::Free, is...) =
    getindex(x, map(f -> f isa Free ? f : Free(f), is)...)
Base.getindex(x, y::Free) = Free(x[y.x], gens(y))

# special case for integer literals
function Base.getindex(x::Free, is::AbstractVector)
    T = eltype(is)
    if x.x isa List && (T <: Integer ||
                        typeintersect(T, Integer) !== Union{} &&
                        all(x -> isa(x, Integer), is))
        Free(List(x.x.xs[is]), gens(x))
    else
        getindex(x, Free(is))
    end
end


## Lazy

gens(l::Lazy) = sort!(unique!(pushgens!(Symbol[], l)::Vector{Symbol}))

Base.:(==)(k::Lazy, l::Lazy) = false

Lazy(x::Lazy) = x
Lazy(x::AbstractSLProgram) = compile(Lazy, x)
Lazy(x) = Const(x)

evaluate(l::Lazy, xs) = evaluate(gens(l), l, xs)
# TODO: remove the 3-arg evaluate methods?


### Const

struct Const{T} <: Lazy
    c::T

    # TODO: this is a hack, delete ASAP
    Const{T}(x::AbstractSLProgram) where {T} = new{T}(x)
    Const(x) = new{typeof(x)}(x)
end

Const(x::AbstractSLProgram) = Const{typeof(x)}(x)


Base.show(io::IO, c::Const) = print(io, c.c)

pushgens!(gs, c::Const) = gs
constantstype(l::Const{T}) where {T} = T

Base.:(==)(k::Const, l::Const) = k.c == l.c

evaluate(gs, c::Const, xs) = c.c

maxinput(c::Const) = 0


### Input

struct Input <: Lazy
    n::Int
end

function Base.show(io::IO, i::Input)
    syms = io[:slp_free_gens]
    print(io, syms[i.n])
end

evaluate(gs, i::Input, xs) = xs[i.n]

maxinput(i::Input) = i.n

Base.:(==)(i::Input, j::Input) = i.n == j.n

constantstype(::Input) = Union{}


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

maxinput(g::Gen) = throw(ArgumentError("logic error: Gen not allowed in Free"))


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

maxinput(p::Plus) = mapreduce(maxinput, max, p.xs)


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

maxinput(m::Minus) = max(maxinput(m.p), maxinput(m.q))


### UniMinus

struct UniMinus <: Lazy
    p::Lazy
end

Base.show(io::IO, p::UniMinus) = print(io, "(-", p.p, ')')

pushgens!(gs, l::UniMinus) = pushgens!(gs, l.p)

constantstype(p::UniMinus) = constantstype(p.p)

Base.:(==)(k::UniMinus, l::UniMinus) = k.p == l.p

evaluate(gs, p::UniMinus, xs) = -evaluate(gs, p.p, xs)

maxinput(p::UniMinus) = maxinput(p.p)


### Times

struct Times <: Lazy
    xs::Vector{Lazy}
end

Times(xs::Lazy...) = Times(collect(Lazy, xs))

function Base.show(io::IO, p::Times)
    print(io, '(')
    join(io, p.xs, '*')
    print(io, ')')
end

pushgens!(gs, p::Times) = foldl(pushgens!, p.xs, init=gs)

constantstype(p::Times) =
    mapreduce(constantstype, typejoin, p.xs, init=Union{})

Base.:(==)(k::Times, l::Times) = k.xs == l.xs

evaluate(gs, p::Times, xs) = mapreduce(q -> evaluate(gs, q, xs), *, p.xs)

maxinput(p::Times) = mapreduce(maxinput, max, p.xs)


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

maxinput(e::Exp) = maxinput(e.p)


### Decision

struct Decision <: Lazy
    ps::Vector{Tuple{Lazy,Int}}
end

test(p::Lazy, n) = Decision(Tuple{Lazy,Int}[(p, n)])

Base.show(io::IO, d::Decision) =
    join(io, (sprint(print, "test(", p, ", ", n, ")", context=io)
              for (p, n) in d.ps),
         " & ")

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

maxinput(p::Decision) = mapreduce(x -> maxinput(x[1]), max, p.ps)


### List

struct List <: Lazy
    xs::Vector{Lazy}
end

function Base.show(io::IO, l::List)
    print(io, "list([")
    join(io, l.xs, ", ")
    print(io, "])")
end

constantstype(l::List) =
    mapreduce(constantstype, typejoin, l.xs, init=Union{})

pushgens!(gs, l::List) = foldl(pushgens!, l.xs, init=gs)

Base.:(==)(p::List, q::List) = p.xs == q.xs

evaluate(gs, l::List, xs) =
    list(eltype(xs)[evaluate(gs, p, xs) for p in l.xs])

maxinput(l::List) = mapreduce(maxinput, max, l.xs)


### Compose

struct Compose <: Lazy
    p::Lazy
    q::Lazy # q must return a list!

    function Compose(p::Lazy, q::Lazy)
        q isa List || throw(ArgumentError(
            "second argument of Compose must return a list"))
        new(p, q)
    end
end

Base.show(io::IO, c::Compose) = print(io, c.p, " ∘ ", c.q)

constantstype(c::Compose) = typejoin(constantstype(c.p, constantstype(c.q)))

pushgens!(gs, c::Compose) = pushgens!(gs, c.q)

Base.:(==)(k::Compose, l::Compose) = k.p == l.p && k.q == l.q

evaluate(gs, p::Compose, xs) = evaluate(gs, p.p, evaluate(gs, p.q, xs))

maxinput(m::Compose) = maxinput(m.q)


### Getindex

struct Getindex <: Lazy
    p::Lazy
    is::Vector{Lazy}
end

function Base.show(io::IO, g::Getindex)
    print(io, g.p, '[')
    join(io, g.is, ", ")
    print(io, ']')
end

pushgens!(gs, l::Getindex) = foldl(pushgens!, l.is, init=pushgens!(gs, l.p))

constantstype(m::Getindex) = mapreduce(constantstype, typejoin, m.is, init=constantstype(m.p))

Base.:(==)(k::Getindex, l::Getindex) = k.p == l.p && k.is == l.is

evaluate(gs, p::Getindex, xs) = evaluate(gs, p.p, xs)[(evaluate(gs, i, xs) for i in p.is)...]

maxinput(m::Getindex) = mapreduce(maxinput, max, m.is, init=maxinput(m.p))


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

#### getindex

Base.getindex(x::Lazy, i::Lazy...) = Getindex(x, collect(i))


### adhoc binary ops

*(x, y::Lazy) = Const(x) * y
*(x::Lazy, y) = x * Const(y)

+(x, y::Lazy) = Const(x) + y
+(x::Lazy, y) = x + Const(y)

-(x, y::Lazy) = Const(x) - y
-(x::Lazy, y) = x - Const(y)

Base.getindex(x, i::Lazy) = Getindex(Const(x), [i])
Base.getindex(x::Lazy, is...) = Getindex(x, Lazy[i isa Lazy ? i : Const(i) for i in is])


## compile to SLProgram

# TODO: this is legacy, only for tests
SLProgram(l::Lazy) = compile(SLProgram, l)
SLProgram{T}(l::Lazy) where {T} = compile(SLProgram{T}, l)

compile(::Type{SLProgram}, l::Lazy) = compile(SLProgram{constantstype(l)}, l)

function compile(::Type{SLProgram{T}}, l::Lazy, gs=gens(l)) where T
    p = SLProgram{T}()
    i = pushlazy!(p, l, gs)
    pushfinalize!(p, i)
end

pushlazy!(p::SLProgram, l::Const, gs) = pushconst!(p, l.c)

pushlazy!(p::SLProgram, l::Input, gs) = input(l.n)

pushlazy!(p::SLProgram, l::Gen, gs) = input(findfirst(==(l.g), gs))

function pushlazy!(p, l::Union{Plus,Times}, gs)
    # TODO: handle isempty(p.xs) ?
    op = l isa Plus ? plus : times
    x, xs = Iterators.peel(l.xs)
    i = pushlazy!(p, x, gs)
    for x in xs
        j = pushlazy!(p, x, gs)
        i = pushop!(p, op, i, j)
    end
    i
end

function pushlazy!(p, l::Minus, gs)
    i = pushlazy!(p, l.p, gs)
    j = pushlazy!(p, l.q, gs)
    pushop!(p, minus, i, j)
end

pushlazy!(p, l::UniMinus, gs) = pushop!(p, uniminus, pushlazy!(p, l.p, gs))

pushlazy!(p, l::Exp, gs) =
    pushop!(p, exponentiate, pushlazy!(p, l.p, gs), intarg(l.e))

function pushlazy!(p, l::Decision, gs)
    local k
    for (x, i) in l.ps
        d = pushlazy!(p, x, gs)
        k = pushop!(p, decision, d, pushint!(p, i))
    end
    k
end

function pushlazy!(p, g::Getindex, gs)
    length(g.is) == 1 ||
        throw(ArgumentError("getindex with multiple indices not supported"))
    i = pushlazy!(p, g.p, gs)
    j = pushlazy!(p, g.is[1], gs)
    pushop!(p, getindex_, i, j)
end
