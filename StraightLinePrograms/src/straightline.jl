## constants

const opmask      = 0xff00000000000000
const argmask     = 0x000000000fffffff
const inputmark   = 0x0000000008000000
const cstmark     = 0x0000000004000000
const intmark     = inputmark | cstmark
const typemark    = intmark
const payloadmask = 0x0000000003ffffff
const negbit      = 0x0000000002000000
const argshift    = 28 # argmask == 2^argshift - 1

@assert payloadmask == argmask ⊻ inputmark ⊻ cstmark

const dectrue = argmask

import Nemo


## Op, Line, Arg

struct Op
    x::UInt64
end

struct Line
    x::UInt64
end

struct Arg
    x::UInt64

    function Arg(x::UInt64)
        iszero(x & ~argmask ) ||
            throw(ArgumentError("argument too big"))
        new(x)
    end
end

Arg(x::Integer) = Arg(UInt64(x))


## predicates

const showop = Dict{Op,String}()

for (i, (op, unary, showchar)) in enumerate([(:assign       , false, "->"),
                                             (:uniminus     , true,  "-"),
                                             (:plus         , false, "+"),
                                             (:minus        , false, "-"),
                                             (:times        , false, "*"),
                                             (:divide       , false, "/"),
                                             (:exponentiate , true,  "^"),
                                             (:keep         , true,  "keep"),
                                             (:decision     , false, "&"),
                                             (:getindex_    , false, "[]"),
                                             ])
    isop = Symbol(:is, op)
    c = UInt64(i) << (2*argshift)
    if unary
        c |= 0x8000000000000000
    end
    @eval begin
        const $op = Op($c)
        $isop(op::Op) = op === $op
    end
    showop[Op(c)] = showchar
end

isquasiunary(op) = (op.x & 0x8000000000000000) != 0
isunary(op) = isquasiunary(op) & (op != exponentiate)


## raw manips

Line(op::Op, i::Arg, j::Arg) = Line(op.x | i.x << argshift | j.x)

pack(op::Op, i::Arg, j::Arg) = Line(op, i, j)

function unpack(line::Line)
    line = line.x
    op = opmask & line
    j = line & argmask
    i = (line >> argshift) & argmask
    Op(op), Arg(i), Arg(j)
end


## show

function Base.show(io::IO, l::Line)
    op, i, j = unpack(l)
    print(io, op, " :  ", i, " , ", j)
end

Base.show(io::IO, op::Op) = print(io, showop[op])

function Base.show(io::IO, x::Arg)
    if isint(x)
        print(io, 'i', intidx(x))
    elseif isinput(x)
        print(io, '$', inputidx(x))
    elseif isconstant(x)
        print(io, '+', constantidx(x))
    else
        print(io, ' ', x.x)
    end
end


## SLProgram

mutable struct SLProgram{T} <: AbstractSLProgram
    cs::Vector{T}       # constants
    lines::Vector{Line} # instructions
    int::Int            # number of stored Int at the beginning of lines
    len::Int            # length of result vector
                        # (where to store next result for any new line, minus 1)
    ret::Arg            # result index to return (0 for all)
    f::Ref{Function}    # compiled execution
end

SLProgram{T}() where {T} = SLProgram(T[], Line[], 0, 0, Arg(0), Ref{Function}())
SLProgram() = SLProgram{Union{}}()

# return an input
function SLProgram{T}(i::Integer) where {T}
    p = SLProgram{T}()
    pushfinalize!(p, input(i))
end

SLProgram(i::Integer) = SLProgram{Union{}}(i)

function Base.copy!(q::SLProgram, p::SLProgram)
    copy!(q.cs, p.cs)
    copy!(q.lines, p.lines)
    q.int = p.int
    q.len = p.len
    q.ret = p.ret
    q
end

copy_oftype(p::SLProgram, ::Type{T}) where {T} = copy!(SLProgram{T}(), p)

Base.copy(p::SLProgram{T}) where {T} = copy_oftype(p, T)

Base.:(==)(p::SLProgram, q::SLProgram) =
    p.cs == q.cs && p.lines == q.lines && p.int == q.int && p.ret == q.ret

constants(p::SLProgram) = p.cs
_integers(p::SLProgram) = @view p.lines[1:p.int]
integers(p::SLProgram) = @view p.lines[p.int:-1:1]
lines(p::SLProgram) = @view p.lines[1+p.int : end]
linesindices(p::SLProgram) = 1+p.int:lastindex(p.lines)
nsteps(p::SLProgram) = length(p.lines) - p.int

constantstype(p::SLProgram{T}) where {T} = T

# return the (max) number of inputs
function ninputs(p::SLProgram)
    m0 = isinput(p.ret) ? inputidx(p.ret) % Int : 0
    m1 = mapreduce(max, lines(p); init=0) do line
        op, i, j = unpack(line)
        max(isinput(i) ? inputidx(i) : 0,
            isinput(j) ? inputidx(j) : 0) % Int
    end
    max(m0, m1)
end

slpgen(n::Integer) = SLProgram(n)
slpgens(n::Integer) = [SLProgram(i) for i=1:n]

gens(::Type{SLProgram}, n::Integer) = slpgens(n)
gens(::Type{SLProgram{T}}, n::Integer) where {T} = [SLProgram{T}(i) for i=1:n]

slpcst(c) = SLProgram(Const(c))


## show

# old basic version, useful when normal show is broken
function showsimple(io::IO, p::SLProgram)
    println("SLProgram:")
    if !isempty(constants(p))
        println("with constants:")
        for (i, c) in enumerate(constants(p))
            println(io, i, " | ", c)
        end
    end
    if !isempty(integers(p))
        println("with integers:")
        for (i, c) in enumerate(integers(p))
            println(io, i, " | ", c.x % Int)
        end
    end
    if !isempty(lines(p))
        println("with lines:")
        for (i, c) in enumerate(lines(p))
            println(io, i, " | ", c)
        end
    end
    println("return: ", p.ret == Arg(dectrue) ? "true" : p.ret)
end

showsimple(p::SLProgram) = showsimple(stdout, p)

slpsyms(n::Integer) =
    n <= 3 ?
        (:x, :y, :z)[1:n] :
        (Symbol(:x, i) for i = 1:n)

function Base.show(io::IO, p::SLProgram)
    gs = get(io, :SLPsymbols, slpsyms(ninputs(p)))
    show(io, evaluate(p, lazyrecgens(gs), Const))
end

function Base.show(io::IO, ::MIME"text/plain", p::SLProgram{T}) where T
    n = length(lines(p))
    gs = get(io, :SLPsymbols, slpsyms(ninputs(p)))
    syms = lazyrecgens(gs)
    reslazy = Any[]
    if n == 0 && !hasmultireturn(p)
        # trivial program, show only result
        return show(io, retrieve(integers(p), constants(p), syms, reslazy, p.ret))
    end

    str(x...) = sprint(print, x...; context=io)
    strlines = Vector{String}[]
    widths = [0, 0, 0, 0, 0]

    ptmp = SLProgram{T}()
    copy!(ptmp.cs, p.cs)
    copy!(ptmp.lines, _integers(p))
    ptmp.int = p.int

    cs = constants(p)
    ints = integers(p)

    showarg(x) = string(retrieve(ints, cs, syms,
                                 ['#'*string(i) for i in eachindex(lines(p))],
                                 x))

    for line in lines(p)
        op, i, j = unpack(line)
        k = updatelen!(ptmp, op, i, j)
        push!(ptmp.lines, line)
        pushfinalize!(ptmp, Arg(k))

        line = String[]
        push!(strlines, line)
        if iskeep(op)
            push!(line, "keep:")
            if i.x == 0
                push!(line, " nothing")
            elseif i.x == 1
                push!(line, " #1")
            else
                push!(line, " #1..#$(i.x)")
            end
            continue
        end

        x = showarg(i)
        y = isunary(op) || isassign(op) ? "" :
            isexponentiate(op) ? string(getint(j)) :
            isquasiunary(op) ? string(j.x) :
            showarg(j)

        if isdecision(op)
            push!(line, "test:")
            push!(line, " order($x) == $y || return false")
            continue
        end

        push!(line, str('#', string(k), " ="))
        # TODO: add a way to evaluate step by step instead of re-evaluating
        # from the beginning each time
        plazy = evaluate!(reslazy, ptmp, syms, Const)

        strop = isassign(op) ? "" : showop[op]
        push!(line, str(strop), x, y, str(plazy))
        widths .= max.(widths, textwidth.(line))
    end
    for (k, line) in enumerate(strlines)
        k == 1 || println(io)
        if line[1] ∈ ["keep:", "test:"]
            join(io, line)
            continue
        end
        print(io, ' '^(widths[1]-length(line[1])), line[1])
        for i = 2:4
            print(io, ' '^(1+widths[i]-textwidth(line[i])), line[i])
        end
        print(io, "  ==>  ", line[5])
    end

    print(io, "\nreturn: ")
    # TODO: ^--- remove when unnecessary?
    if hasmultireturn(p)
        print(io, '[')
        join(io,
             (showarg(Arg(i)) for i in 1:p.len),
             ", ")
        print(io, ']')
    elseif hasdecision(p)
        print(io, true)
    else
        print(io, showarg(p.ret))
    end

end


## building SLProgram

isregister(i::Arg) = typemark & i.x == 0
registeridx(i::Arg) = i.x
asregister(i::Integer) = Arg(UInt64(i)) # TODO: sanitize input

# return #ref for i-th input
function input(i::Integer)
    i = Int(i)
    @assert i < cstmark
    Arg((i % UInt64) | inputmark)
end

isinput(i::Arg) = typemark & i.x === inputmark
inputidx(i::Arg) = i.x ⊻ inputmark

isconstant(i::Arg) = typemark & i.x === cstmark
constantidx(i::Arg) = i.x ⊻ cstmark

asconstant(i::Integer) = Arg(UInt64(i) | cstmark)

isint(i::Arg) = typemark & i.x === intmark
asint(i::Int) = Arg(i % UInt64 | intmark)
intidx(i::Arg) = i.x ⊻ intmark

hasmultireturn(p::SLProgram) = p.ret.x == 0
setmultireturn!(p::SLProgram) = (p.ret = Arg(0); p)

hasdecision(p::SLProgram) = p.ret.x == argmask
setdecision!(p::SLProgram) = (p.ret = Arg(argmask); p)


# make an Arg out of x, considered as "computational" integer (e.g. not an index)
# call to create an integer exponent
function intarg(x::Integer)
    y = Int(x) % UInt64
    if x >= 0
        y > (payloadmask ⊻ negbit) &&
            throw(ArgumentError("positive integer argument too big"))
        Arg(y)
    else
        y | (payloadmask >> 1) == typemax(UInt64) || # >> 1 for signbit
            throw(ArgumentError("negative integer argument too small"))
        Arg(y & payloadmask)
    end
end

# retrieve a "computational" integer from an Arg
function getint(x::Arg)
    y = x.x
    iszero(y & ~payloadmask) ||
        throw(ArgumentError("Arg does not contain an Int"))
    if iszero(y & negbit)
        y % Int
    else
        (y | ~payloadmask) % Int
    end
end

# call before mutating, unless p is empty (opposite of pushfinalize!)

# TODO: pushinit! & pushfinalize! not really nececessary anymore?
pushinit!(p::SLProgram) = p.ret

function pushconst!(p::SLProgram{T}, c::T) where T
    push!(constants(p), c)
    l = lastindex(constants(p))
    @assert l < cstmark
    asconstant(l)
end

function pushint!(p::SLProgram, i::Integer)
    pushfirst!(p.lines, Line(Int(i) % UInt64))
    p.int += 1
    @assert p.int < cstmark
    asint(p.int)
end

function updatelen!(p, op, i, j)
    if isassign(op) && j != Arg(0)
        ptr = Int(j.x)
        if ptr == p.len + 1
            p.len += 1
        end
        1 <= ptr <= p.len ||
            throw(ArgumentError("invalid `assign` destination"))
    elseif iskeep(op)
        ptr = Int(i.x)
        ptr <= p.len || throw(ArgumentError("cannot `keep` so many items"))
        p.len = ptr
    elseif isdecision(op)
        ptr = argmask # cf. hasdecision
    else
        ptr = p.len += 1
        @assert ptr < cstmark
    end
    ptr
end

function pushop!(p::SLProgram, op::Op, i::Arg, j::Arg=Arg(0))
    @assert i.x <= argmask && j.x <= argmask
    ptr = updatelen!(p, op, i, j)
    push!(p.lines, Line(op, i, j))
    Arg(ptr % UInt64)
end

function pushfinalize!(p::SLProgram, ret::Arg)
    p.ret = ret
    p
end

function _combine!(p::SLProgram, q::SLProgram; makeinput=identity)
    i1 = pushinit!(p)
    i2 = pushinit!(q)
    koffset = length(constants(p))
    ioffset = p.int
    len = length(lines(p))
    loffset = len - count(lines(p)) do line
        op, i, j = unpack(line)
        isdecision(op)
        # TODO: probably needs also to account for "keep" lines
    end
    p.int += q.int # must be *after* computing len
    prepend!(p.lines, _integers(q))
    append!(p.lines, lines(q))
    append!(constants(p), constants(q))

    @assert length(constants(p)) < cstmark # TODO: should not be @assert

    for n = len+1:lastindex(lines(p))
        op, i, j = unpack(lines(p)[n])
        if isconstant(i)
            i = Arg(i.x + koffset)
        elseif isinput(i)
            i = makeinput(i)
        elseif isint(i)
            i = Arg(i.x + ioffset)
        else
            i = Arg(i.x + loffset)
        end
        if isconstant(j)
            j = Arg(j.x + koffset)
        elseif isinput(j)
            j = makeinput(j)
        elseif isint(j)
            j = Arg(j.x + ioffset)
        elseif !isquasiunary(op)
            # TODO: normalize assign so that j.x is never 0
            if !(isassign(op) && j.x == 0)
                j = Arg(j.x + loffset)
            end
        end
        lines(p)[n] = Line(op, i, j)
        # TODO: write conditionally only when modifications
    end
    if isconstant(i2)
        i2 = Arg(i2.x + koffset)
    elseif isinput(i2)
        i2 = makeinput(i2)
    elseif hasdecision(q)
    elseif hasmultireturn(q)
    elseif isint(i2)
        i2 = Arg(i2.x + ioffset)
    else
        i2 = Arg(i2.x + loffset)
    end
    p.len += q.len
    i1, i2
end

function combine!(op::Op, p::SLProgram, q::SLProgram)
    i = pushop!(p, op, _combine!(p, q)...)
    pushfinalize!(p, i)
end

function combine!(op::Op, p::SLProgram)
    i = pushinit!(p)
    i = pushop!(p, op, i)
    pushfinalize!(p, i)
end

function combine!(op::Op, p::SLProgram, e::Integer)
    i = pushinit!(p)
    i = pushop!(p, op, i, intarg(e))
    pushfinalize!(p, i)
end

applyperm(perm, i, usegetindex) = usegetindex ? perm[i] : perm(i)

function permute_inputs!(p::SLProgram, perm, usegetindex=true)
    for l in linesindices(p)
        op, i, j = unpack(p.lines[l])
        changed = false
        if isinput(i)
            i = input(applyperm(perm, inputidx(i) % Int, usegetindex))
            changed = true
        end
        if isinput(j)
            j = input(applyperm(perm, inputidx(j) % Int, usegetindex))
            changed = true
        end
        if changed
            p.lines[l] = pack(op, i, j)
        end
    end
    p
end


### adhoc

function combine(op::Op, p::SLProgram, q)
    r = copy_oftype(p, typejoin(constantstype(p), typeof(q)))
    i = pushinit!(r)
    j = pushop!(r, op, i, pushconst!(r, q))
    pushfinalize!(r, j)
end

function combine(op::Op, p, q::SLProgram)
    r = copy_oftype(q, typejoin(constantstype(q), typeof(p)))
    i = pushinit!(r)
    j = pushop!(r, op, pushconst!(r, p), i)
    pushfinalize!(r, j)
end


## mutating ops

addeq!(p::SLProgram, q::SLProgram) = combine!(plus, p, q)

subeq!(p::SLProgram, q::SLProgram) = combine!(minus, p, q)

function subeq!(p::SLProgram)
    combine!(uniminus, p)
    p
end

muleq!(p::SLProgram, q::SLProgram) = combine!(times, p, q)

function expeq!(p::SLProgram, e::Integer)
    combine!(exponentiate, p, e)
    p
end

function testeq!(p::SLProgram, q::SLProgram)
    hasdecision(p) && hasdecision(q) || throw(ArgumentError(
        "cannot &-combine two programs which are not decisions"))
    _combine!(p, q)
    setdecision!(p)
end


## unary/binary ops

copy_jointype(p::SLProgram, q::SLProgram) =
    copy_oftype(p, typejoin(constantstype(p), constantstype(q)))

+(p::SLProgram, q::SLProgram) = addeq!(copy_jointype(p, q), q)

*(p::SLProgram, q::SLProgram) = muleq!(copy_jointype(p, q), q)

-(p::SLProgram, q::SLProgram) = subeq!(copy_jointype(p, q), q)

-(p::SLProgram) = subeq!(copy(p))

^(p::SLProgram, e::Integer) = expeq!(copy(p), e)

Base.literal_pow(::typeof(^), p::SLProgram, ::Val{e}) where {e} = p^e

Base.:&(p::SLProgram, q::SLProgram) = testeq!(copy_jointype(p, q), q)

function test!(p::SLProgram, x::Integer)
    hasdecision(p) && throw(ArgumentError("SLProgram is already a decision"))
    i = pushinit!(p)
    j = pushint!(p, x)
    pushfinalize!(p, pushop!(p, decision, i, j))
end

test(p::SLProgram, x::Integer) = test!(copy(p), x)

function list(::Type{SL}, ps) where {SL<:SLProgram}
    idx = Arg[]
    q = SL === SLProgram ? SLProgram{Any}() : SL()
    # TODO: ^^^--- use typejoin or smthg instead of Any
    for (n, p) in enumerate(ps)
        _, i = _combine!(q, p)
        if !isregister(i)
            # we write explicitly in res-vector, so that the indices increase
            # strictly and it's always safe in the final loop to copy from
            # this `i` into final destination
            i = pushop!(q, assign, i)
            @assert isregister(i)
        end
        push!(idx, i)
        @assert nsteps(q) >= registeridx(i)
    end
    for (n, i) in enumerate(idx)
        if registeridx(i) > n
            pushop!(q, assign, i, Arg(n))
        else
            @assert registeridx(i) == n
        end
    end
    if nsteps(q) > length(idx)
        pushop!(q, keep, Arg(length(idx)))
    end
    q
end

function composewith!(q::SLProgram, p::SLProgram)
    hasmultireturn(q) || throw(ArgumentError(
        "second argument of `compose` must return a list"))
    ninp = q.len # max ninputs that p can use
    _, j = _combine!(q, p, makeinput = function (i::Arg)
                  idx = inputidx(i)
                  idx <= ninp || throw(ArgumentError(
                      "inner compose argument does not provide enough input"))
                  asregister(idx)
              end)
    if j.x == 0 && p.len != 0 # multireturn
        newlen = q.len - ninp
        for i = 1:newlen
            pushop!(q, assign, asregister(i+ninp), asregister(i))
        end
        pushop!(q, keep, asregister(newlen))
    end
    pushfinalize!(q, j)
end


compose(p::SLProgram, q::SLProgram) = composewith!(copy_jointype(q, p), p)

getindex!(p::SLProgram, q::SLProgram) = combine!(getindex_, p, q)
Base.getindex(p::SLProgram, q::SLProgram) = getindex!(copy_jointype(p, q), q)


### adhoc

+(p::SLProgram, q) = combine(plus, p, q)
+(p, q::SLProgram) = combine(plus, p, q)
*(p::SLProgram, q) = combine(times, p, q)
*(p, q::SLProgram) = combine(times, p, q)
-(p::SLProgram, q) = combine(minus, p, q)
-(p, q::SLProgram) = combine(minus, p, q)

function getindex!(p::SLProgram, i::Integer)
    if hasmultireturn(p)
        pushfinalize!(p, asregister(i))
    else
        a = pushinit!(p)
        k = pushop!(p, getindex_, a, pushint!(p, i))
        pushfinalize!(p, k)
    end
end

Base.getindex(p::SLProgram, i::Integer) = getindex!(copy(p), i)

## conversion SLProgram -> LazyRec

lazyrecgens(gs) = Any[Gen(s) for s in gs]
lazyrecgens(n::Integer) = lazyrecgens(slpsyms(n))

# strictly convenience function
aslazyrec(p::SLProgram, gs=lazyrecgens(ninputs(p))) = LazyRec(evaluate(p, gs))


## evaluate

getres(p::SLProgram, xs) = Vector{eltype(xs)}(undef, length(p.lines))
# length(p.lines) might be an overestimate in some cases, but
# this generally limits the number of allocations

function evaluate(p::SLProgram{T}, xs::Vector{S}, conv::F=identity
                  ) where {T,S,F}
    if isassigned(p.f)
      #not too good, but avoids the worldage problem and, evetually the overhead of
      #invokelatest
      try
        return p.f[](xs)::S
      catch 
        return Base.invokelatest(p.f[], xs)::S
      end
    else
        evaluate!(getres(p, xs), p, xs, conv)
    end
end

function evaluates(p::SLProgram{T}, xs::Vector{S}, conv::F=identity
                   ) where {T,S,F}
    res = getres(p, xs)
    evaluate!(res, p, xs, conv)
    res
end

retrieve(ints, cs, xs, res, i, conv::F=identity) where {F} =
    isint(i)      ? ints[intidx(i)].x % Int :
    isconstant(i) ? conv(cs[constantidx(i)]) :
    isinput(i)    ? xs[inputidx(i)] :
    res[i.x]

function evaluate!(res::Vector{S}, p::SLProgram{T}, xs::Vector{S},
                   conv::F=identity) where {S,T,F}
    # TODO: handle isempty(lines(p))
    # TODO: use inplace (addeq!, mul!, ... ) when applicable
    # TODO: add permutation of input?
    empty!(res)

    cs = constants(p)
    ints = integers(p)

    local decide::Union{S,Bool}

    for line in lines(p)
        local r::S
        op, i, j = unpack(line)
        if iskeep(op)
            @assert isregister(i)
            resize!(res, i.x)
            continue
        end
        x = retrieve(ints, cs, xs, res, i, conv)
        if isexponentiate(op)
            r = x^getint(j) # TODO: support bigger j
        elseif isassign(op)
            dst = j.x % Int
            if dst != 0 && dst != lastindex(res) + 1
                res[dst] = x
                continue
            end
            r = x
        elseif isuniminus(op)
            r = -x
        else
            y = retrieve(ints, cs, xs, res, j, conv)
            if isplus(op)
                r = x + y
            elseif isminus(op)
                r = x - y
            elseif istimes(op)
                r = x * y
            elseif isdivide(op)
                r = divexact(x, y)
            elseif isgetindex_(op)
                r = x[y]
            elseif isdecision(op)
                t = test(x, y)
                t === false && return false
                if @isdefined(decide)
                    decide &= t
                else
                    decide = t
                end
                continue
            else
                throw(ArgumentError("unknown operation"))
            end
        end
        push!(res, r)
    end

    @assert length(res) == p.len
    if hasdecision(p)
        decide
    elseif hasmultireturn(p)
        list(res)
    else
        retrieve(ints, cs, xs, res, p.ret, conv)
    end
end

function order end
# to be specialized by users

test(x, o) = order(x) == o

list(res) = list(eltype(res), res)

list(::Type{T}, res) where {T} = res


## compile!

cretrieve(i) =
    isinput(i) ? Symbol(:x, inputidx(i)) => inputidx(i) :
    isconstant(i) ? Symbol(:c, constantidx(i)) => -1 :
    Symbol(:res, i.x) => 0

# TODO: handle the "conv" argument like in evaluate!
# (works without, but there can be type-instability)

# return compiled execution function f, and updates
# p.f[] = f, which is not invalidated when p is mutated
function compile!(p::SLProgram; isPoly::Bool = false)
    res = Expr[]
    fn = :(function (xs::Vector{T}) where {T}
           end)
    k = 0
    cs = constants(p)
    for k in eachindex(cs)
        push!(res, :($(Symbol(:c, k)) = @inbounds $cs[$k]))
    end
    mininput = 0
    for line in lines(p)
        k += 1
        rk = Symbol(:res, k)
        op, i, j = unpack(line)
        x, idx = cretrieve(i)
        mininput = max(mininput, idx)
        line =
            if isexponentiate(op)
                :($rk = $x^$(getint(j)))
            elseif isassign(op)
                :($rk = $x)
            elseif isuniminus(op)
                :($rk = -$x)
            else
                y, idy = cretrieve(j)
                mininput = max(mininput, idy)
                if isplus(op)
                    if idx == 0 && isPoly
                      :($rk = Nemo.addeq!($x, $y))
                    else
                      :($rk = $x + $y)
                    end
                elseif isminus(op)
                    if idx == 0 && isPoly
                      :($rk = Nemo.subeq!($x, $y))
                    else
                      :($rk = $x - $y)
                    end
                elseif istimes(op)
                    if idx == 0 && isPoly
                      :($rk = Nemo.mul!($x, $x, $y))
                    else
                      :($rk = $x * $y)
                    end
                elseif isdivide(op)
                    :($rk = divexact($x, $y))
                end
            end
        push!(res, line)
    end
    for k = 1:mininput-1
        pushfirst!(res, :($(Symbol(:x, k)) = @inbounds xs[$k ]))
    end
    if mininput >= 1
        pushfirst!(res, :($(Symbol(:x, mininput)) = @inbounds xs[$mininput]))
    end
    append!(fn.args[2].args, res)
    global last_bla = fn
    p.f[] = eval(fn)
end
