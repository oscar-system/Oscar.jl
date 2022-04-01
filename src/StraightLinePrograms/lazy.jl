## Lazy

abstract type LazyRec end

struct Lazy <: AbstractSLProgram
    x::LazyRec # must not contain Gen
    gens::Vector{Symbol}
end

Lazy(x::LazyRec) = Lazy(x, collect(Symbol, slpsyms(maxinput(x))))

Lazy(x) = Lazy(Const(x), Symbol[])

function lazygens(n::Integer, syms=slpsyms(n))
    if !isa(syms, Vector{Symbol})
        syms = collect(Symbol, syms)
    end
    Lazy[Lazy(Input(i), syms) for i=1:n]
end

lazygens(syms::AbstractVector{Symbol}) = lazygens(length(syms), syms)

gens(::Type{Lazy}, n::Integer) = lazygens(n)

ngens(f::Lazy) = length(f.gens)

gens(f::Lazy) = f.gens

Base.show(io::IO, f::Lazy) =
    show(IOContext(io, :slp_lazy_gens => gens(f)), f.x)

evaluate(f::Lazy, xs) = evaluate(gens(f), f.x, xs, IdDict{LazyRec,Any}())

# check compatibility and return the biggest of the two gens arrays
function gens(x::Lazy, xs::Lazy...)
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
Base.:(==)(x::Lazy, y::Lazy) = x.x == y.x

compile(::Type{SLProgram}, f::Lazy) =
    compile(SLProgram{constantstype(f.x)}, f)

function compile(::Type{SLProgram{T}}, f::Lazy) where T
    p = SLProgram{T}()
    i = pushlazy!(p, f.x, gens(f))
    pushfinalize!(p, i)
end

compile(f::Lazy) = compile(SLProgram, f)


### unary/binary ops

+(x::Lazy, y::Lazy) = Lazy(x.x + y.x, gens(x, y))
-(x::Lazy, y::Lazy) = Lazy(x.x - y.x, gens(x, y))
*(x::Lazy, y::Lazy) = Lazy(x.x * y.x, gens(x, y))

Base.:(&)(x::Lazy, y::Lazy) = Lazy(x.x & y.x, gens(x, y))

-(x::Lazy) = Lazy(-x.x, gens(x))
^(x::Lazy, e::Integer) = Lazy(x.x^e, gens(x))

Base.literal_pow(::typeof(^), p::Lazy, ::Val{e}) where {e} = p^e

test(x::Lazy, n) = Lazy(test(x.x, n), gens(x))

# TODO: don't splat
list(::Type{Lazy}, xs) = Lazy(List(LazyRec[x.x for x in xs]), gens(xs...))

function compose(p::Lazy, q::Lazy; flatten=true)
    if flatten
        q.x isa List || throw(ArgumentError(
            "first argument must return a list"))
        gs = gens(q)
        evaluate(p, [Lazy(qi, gs) for qi in q.x.xs])::Lazy
    else
        Lazy(Compose(p.x, q.x), gens(q))
    end
end

Base.getindex(p::Lazy, is::Lazy...) =
    Lazy(Getindex(p.x, LazyRec[i.x for i in is]), gens(p, is...))

call(f, xs...) = Lazy(Call(f, [LazyRec(x) for x in xs]))


#### adhoc

+(x::Lazy, y) = Lazy(x.x + y, gens(x))
+(x, y::Lazy) = Lazy(x + y.x, gens(y))

-(x::Lazy, y) = Lazy(x.x - y, gens(x))
-(x, y::Lazy) = Lazy(x - y.x, gens(y))

*(x::Lazy, y) = Lazy(x.x * y, gens(x))
*(x, y::Lazy) = Lazy(x * y.x, gens(y))

Base.getindex(x::Lazy, is...) =
    getindex(x, map(f -> f isa Lazy ? f : Lazy(f), is)...)
Base.getindex(x, y::Lazy) = Lazy(x[y.x], gens(y))

# special case for integer literals

function Base.getindex(x::Lazy, i::Integer)
    if x.x isa List
        Lazy(x.x.xs[i], gens(x))
    else
        getindex(x, Lazy(i))
    end
end

function Base.getindex(x::Lazy, is::AbstractVector)
    T = eltype(is)
    if x.x isa List && (T <: Integer ||
                        typeintersect(T, Integer) !== Union{} &&
                        all(x -> isa(x, Integer), is))
        Lazy(List(x.x.xs[is]), gens(x))
    else
        getindex(x, Lazy(is))
    end
end


## LazyRec

gens(l::LazyRec) = sort!(unique!(pushgens!(Symbol[], l)::Vector{Symbol}))

Base.:(==)(k::LazyRec, l::LazyRec) = false

LazyRec(x::LazyRec) = x
LazyRec(x::AbstractSLProgram) = compile(LazyRec, x)
LazyRec(x) = Const(x)
LazyRec(x::Lazy) = x.x

evaluate(l::LazyRec, xs) = evaluate(gens(l), l, xs, IdDict{LazyRec,Any}())
# TODO: remove the 3-arg evaluate methods?

const sentinel = 0x54f765833cf932f72a3dd18e0dc1d839

function evaluate(gs, l::LazyRec, xs, dict)
    r = get(dict, l, sentinel)
    r !== sentinel && return r
    dict[l] = _evaluate(gs, l, xs, dict)
end

### Const

struct Const{T} <: LazyRec
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

evaluate(gs, c::Const, xs, dict) = c.c

maxinput(c::Const) = 0


### Input

struct Input <: LazyRec
    n::Int
end

function Base.show(io::IO, i::Input)
    syms = io[:slp_lazy_gens]
    print(io, syms[i.n])
end

evaluate(gs, i::Input, xs, dict) = xs[i.n]

maxinput(i::Input) = i.n

Base.:(==)(i::Input, j::Input) = i.n == j.n

constantstype(::Input) = Union{}


### Gen

struct Gen <: LazyRec
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

# TODO: test if dict should be used (performance only)
evaluate(gs, g::Gen, xs, dict) = xs[findfirst(==(g.g), gs)]

maxinput(g::Gen) = throw(ArgumentError("logic error: Gen not allowed in Lazy"))


### Plus

struct Plus <: LazyRec
    xs::Vector{LazyRec}
end

Plus(xs::LazyRec...) = Plus(collect(LazyRec, xs))

function Base.show(io::IO, p::Plus)
    print(io, '(')
    join(io, p.xs, " + ")
    print(io, ')')
end

pushgens!(gs, p::Plus) = foldl(pushgens!, p.xs, init=gs)

constantstype(p::Plus) =
    mapreduce(constantstype, typejoin, p.xs, init=Union{})

Base.:(==)(k::Plus, l::Plus) = k.xs == l.xs

_evaluate(gs, p::Plus, xs, dict) =
    mapreduce(q -> evaluate(gs, q, xs, dict), +, p.xs)

maxinput(p::Plus) = mapreduce(maxinput, max, p.xs)


### Minus

struct Minus <: LazyRec
    p::LazyRec
    q::LazyRec
end

Base.show(io::IO, p::Minus) = print(io, '(', p.p, " - ", p.q, ')')

pushgens!(gs, l::Minus) = pushgens!(pushgens!(gs, l.p), l.q)

constantstype(m::Minus) = typejoin(constantstype(m.p), constantstype(m.q))

Base.:(==)(k::Minus, l::Minus) = k.p == l.p && k.q == l.q

_evaluate(gs, p::Minus, xs, dict) =
    evaluate(gs, p.p, xs, dict) - evaluate(gs, p.q, xs, dict)

maxinput(m::Minus) = max(maxinput(m.p), maxinput(m.q))


### UniMinus

struct UniMinus <: LazyRec
    p::LazyRec
end

Base.show(io::IO, p::UniMinus) = print(io, "(-", p.p, ')')

pushgens!(gs, l::UniMinus) = pushgens!(gs, l.p)

constantstype(p::UniMinus) = constantstype(p.p)

Base.:(==)(k::UniMinus, l::UniMinus) = k.p == l.p

_evaluate(gs, p::UniMinus, xs, dict) = -evaluate(gs, p.p, xs, dict)

maxinput(p::UniMinus) = maxinput(p.p)


### Times

struct Times <: LazyRec
    xs::Vector{LazyRec}
end

Times(xs::LazyRec...) = Times(collect(LazyRec, xs))

function Base.show(io::IO, p::Times)
    print(io, '(')
    join(io, p.xs, '*')
    print(io, ')')
end

pushgens!(gs, p::Times) = foldl(pushgens!, p.xs, init=gs)

constantstype(p::Times) =
    mapreduce(constantstype, typejoin, p.xs, init=Union{})

Base.:(==)(k::Times, l::Times) = k.xs == l.xs

_evaluate(gs, p::Times, xs, dict) =
    mapreduce(q -> evaluate(gs, q, xs, dict), *, p.xs)

maxinput(p::Times) = mapreduce(maxinput, max, p.xs)


### Exp

struct Exp <: LazyRec
    p::LazyRec
    e::Int
end

Base.show(io::IO, p::Exp) = print(io, p.p, '^', p.e)

pushgens!(gs, l::Exp) = pushgens!(gs, l.p)

constantstype(p::Exp) = constantstype(p.p)

Base.:(==)(k::Exp, l::Exp) = k.p == l.p && k.e == l.e

Base.literal_pow(::typeof(^), p::LazyRec, ::Val{e}) where {e} = Exp(p, e)

_evaluate(gs, p::Exp, xs, dict) = evaluate(gs, p.p, xs, dict)^p.e

maxinput(e::Exp) = maxinput(e.p)


### Decision

struct Decision <: LazyRec
    ps::Vector{Tuple{LazyRec,Int}}
end

test(p::LazyRec, n) = Decision(Tuple{LazyRec,Int}[(p, n)])

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

function _evaluate(gs, p::Decision, xs, dict)
    (d, i), rest = Iterators.peel(p.ps) # p.ps should never be empty
    res = test(evaluate(gs, d, xs, dict), i)
    res === false && return false
    for (d, i) in rest
        r = test(evaluate(gs, d, xs, dict), i)
        r === false && return false
        res &= r
    end
    res
end

maxinput(p::Decision) = mapreduce(x -> maxinput(x[1]), max, p.ps)


### List

struct List <: LazyRec
    xs::Vector{LazyRec}
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

_evaluate(gs, l::List, xs, dict) =
    list(eltype(xs)[evaluate(gs, p, xs, dict) for p in l.xs])

maxinput(l::List) = mapreduce(maxinput, max, l.xs)


### Compose

struct Compose <: LazyRec
    p::LazyRec
    q::LazyRec # q must return a list!

    function Compose(p::LazyRec, q::LazyRec)
        q isa List || throw(ArgumentError(
            "second argument of Compose must return a list"))
        new(p, q)
    end
end

Base.show(io::IO, c::Compose) = print(io, c.p, " ∘ ", c.q)

constantstype(c::Compose) = typejoin(constantstype(c.p, constantstype(c.q)))

pushgens!(gs, c::Compose) = pushgens!(gs, c.q)

Base.:(==)(k::Compose, l::Compose) = k.p == l.p && k.q == l.q

_evaluate(gs, p::Compose, xs, dict) =
    evaluate(gs, p.p, evaluate(gs, p.q, xs, dict), dict)

maxinput(m::Compose) = maxinput(m.q)


### Getindex

struct Getindex <: LazyRec
    p::LazyRec
    is::Vector{LazyRec}
end

function Base.show(io::IO, g::Getindex)
    print(io, g.p, '[')
    join(io, g.is, ", ")
    print(io, ']')
end

pushgens!(gs, l::Getindex) = foldl(pushgens!, l.is, init=pushgens!(gs, l.p))

constantstype(m::Getindex) = mapreduce(constantstype, typejoin, m.is, init=constantstype(m.p))

Base.:(==)(k::Getindex, l::Getindex) = k.p == l.p && k.is == l.is

_evaluate(gs, p::Getindex, xs, dict) =
    evaluate(gs, p.p, xs, dict)[(evaluate(gs, i, xs, dict) for i in p.is)...]

maxinput(m::Getindex) = mapreduce(maxinput, max, m.is, init=maxinput(m.p))


### Call

struct Call <: LazyRec
    f::Any # callable
    args::Vector{LazyRec}
end


function Base.show(io::IO, f::Call)
    print(io, f.f, '(')
    join(io, f.args, ", ")
    print(io, ')')
end

pushgens!(gs, f::Call) = foldl(pushgens!, f.args, init=Symbol[])

constantstype(f::Call) = mapreduce(constantstype, typejoin, f.args, init=Union{})

Base.:(==)(f::Call, g::Call) = f.f == g.f && f.args == g.args

_evaluate(gs, f::Call, xs, dict) =
    f.f(map(t -> evaluate(gs, t, xs, dict), f.args)...)

maxinput(f::Call) = mapreduce(maxinput, max, f.args)


### binary ops

#### +

+(x::LazyRec, y::LazyRec) = Plus(x, y)

function +(x::Plus, y::LazyRec)
    p = Plus(copy(x.xs))
    push!(p.xs, y)
    p
end

function +(x::LazyRec, y::Plus)
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

-(p::LazyRec, q::LazyRec) = Minus(p, q)
-(p::LazyRec) = UniMinus(p)


#### *

*(x::LazyRec, y::LazyRec) = Times(x, y)

function *(x::Times, y::LazyRec)
    p = Times(copy(x.xs))
    push!(p.xs, y)
    p
end

function *(x::LazyRec, y::Times)
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

^(x::LazyRec, e::Integer) = Exp(x, Int(e))

#### getindex

Base.getindex(x::LazyRec, i::LazyRec...) = Getindex(x, collect(i))


### adhoc binary ops

*(x, y::LazyRec) = Const(x) * y
*(x::LazyRec, y) = x * Const(y)

+(x, y::LazyRec) = Const(x) + y
+(x::LazyRec, y) = x + Const(y)

-(x, y::LazyRec) = Const(x) - y
-(x::LazyRec, y) = x - Const(y)

Base.getindex(x, i::LazyRec) = Getindex(Const(x), [i])
Base.getindex(x::LazyRec, is...) = Getindex(x, LazyRec[i isa LazyRec ? i : Const(i) for i in is])


## compile to SLProgram

# TODO: this is legacy, only for tests
SLProgram(l::LazyRec) = compile(SLProgram, l)
SLProgram{T}(l::LazyRec) where {T} = compile(SLProgram{T}, l)

compile(::Type{SLProgram}, l::LazyRec) = compile(SLProgram{constantstype(l)}, l)

function compile(::Type{SLProgram{T}}, l::LazyRec, gs=gens(l)) where T
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
