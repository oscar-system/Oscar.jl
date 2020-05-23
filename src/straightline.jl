## Op, Line, Arg

struct Op
    x::UInt64
end

struct Line
    x::UInt64
end

struct Arg
    x::UInt64
end


## constants & predicates

const opmask    = 0xff00000000000000
const argmask   = 0x000000000fffffff
const inputmark = 0x0000000008000000
const cstmark   = 0x0000000004000000
const argshift  = 28 # argmask == 2^argshift - 1


const showop = Dict{Op,Char}()

for (i, (op, unary, showchar)) in enumerate([(:uniplus      , true,  '+'),
                                             (:uniminus     , true,  '-'),
                                             (:plus         , false, '+'),
                                             (:minus        , false, '-'),
                                             (:times        , false, '*'),
                                             (:divide       , false, '/'),
                                             (:exponentiate , true,  '^')])
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

pack(op::Op, i, j) = Line(op, Arg(i), Arg(j))

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
    if isinput(x)
        print(io, '$', inputidx(x))
    elseif isconstant(x)
        print(io, '+', constantidx(x))
    else
        print(io, ' ', x.x)
    end
end


## SLProgram

struct SLProgram{T}
    cs::Vector{T}       # constants
    lines::Vector{Line} # instructions
    f::Ref{Function}    # compiled execution
end

SLProgram(cs, lines) = SLProgram(cs, lines, Ref{Function}())

SLProgram{T}() where {T} = SLProgram(T[], Line[])

# return an input
function SLProgram{T}(i::Integer) where {T}
    p = SLProgram{T}()
    pushfinalize!(p, pushop!(p, uniplus, input(i)))
end

function SLProgram(c::Const{T}) where {T}
    p = SLProgram{T}()
    pushfinalize!(p, pushop!(p, uniplus, pushconst!(p, c.c)))
end

constants(p::SLProgram) = p.cs
lines(p::SLProgram) = p.lines


## show

function Base.show(io::IO, ::MIME"text/plain", p::SLProgram)
    println("SLProgram with constants:")
    for (i, c) in enumerate(constants(p))
        println(io, i, " | ", c)
    end
    println("and with lines:")
    for (i, l) in enumerate(lines(p))
        i == 1 || println(io)
        print(io, i, " | ", l)
    end
end


## building SLProgram

# return #ref for i-th input
function input(i::Integer)
    i = Int(i)
    @assert i < cstmark
    Arg((i % UInt64) | inputmark)
end

isinput(i::Arg) = inputmark & i.x != 0
inputidx(i::Arg) = i.x ⊻ inputmark

isconstant(i::Arg) = cstmark & i.x != 0
constantidx(i::Arg) = i.x ⊻ cstmark

asconstant(i::Integer) = Arg(UInt64(i) | cstmark)

# call before mutating, unless p is empty (opposite of pushfinalize!)

function pushinit!(p::SLProgram)
    plines = lines(p)
    op, i, j = unpack(plines[end])
    if isuniplus(op) && (isinput(i) || isconstant(i))
        pop!(plines) # discard trivial instruction
        i
    else
        Arg(lastindex(plines) % UInt64)
    end
end

function pushconst!(p::SLProgram{T}, c::T) where T
    push!(constants(p), c)
    l = lastindex(constants(p))
    @assert l < cstmark
    asconstant(l)
end

function pushop!(p::SLProgram, op::Op, i::Arg, j::Arg=Arg(0))
    @assert i.x <= argmask && j.x <= argmask
    push!(lines(p), Line(op, i, j))
    l = lastindex(lines(p))
    @assert l < cstmark
    Arg(l % UInt64)
end

# make state consistent again
function pushfinalize!(p::SLProgram, ret::Arg)
    k = length(constants(p))
    if isinput(ret) || isconstant(ret) || ret.x != lastindex(lines(p))
        # non-trivial return (i.e. no line or not the result of the last line)
        pushop!(p, uniplus, ret)
    end
    p
end

function _combine!(p::SLProgram{T}, q::SLProgram{T}) where T
    i1 = pushinit!(p)
    koffset = length(constants(p))
    len = length(lines(p))
    append!(lines(p), lines(q))
    append!(constants(p), constants(q))

    @assert length(constants(p)) < cstmark # TODO: should not be @assert
    for n = len+1:lastindex(lines(p))
        op, i, j = unpack(lines(p)[n])
        if isconstant(i)
            i = Arg(i.x + koffset)
        elseif isinput(i)
        else
            i = Arg(i.x + len)
        end
        if isconstant(j)
            j = Arg(j.x + koffset)
        elseif isinput(j)
        elseif !isquasiunary(op)
            j = Arg(j.x + len)
        end
        lines(p)[n] = Line(op, i, j)
        # TODO: write conditionally only when modifications
    end
    i2 = pushinit!(p)
    i1, i2
end

function combine!(op::Op, p::SLProgram{T}, q::SLProgram{T}) where {T}
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
    i = pushop!(p, op, i, Arg(UInt64(e)))
    pushfinalize!(p, i)
end


## mutating ops

addeq!(p::SLProgram{T}, q::SLProgram{T}) where {T} = combine!(plus, p, q)

subeq!(p::SLProgram{T}, q::SLProgram{T}) where {T} = combine!(minus, p, q)

function subeq!(p::SLProgram)
    combine!(uniminus, p)
    p
end

muleq!(p::SLProgram{T}, q::SLProgram{T}) where {T} = combine!(times, p, q)

function expeq!(p::SLProgram, e::Integer)
    combine!(exponentiate, p, e)
    p
end


## execute

function execute(p::SLProgram{T}, xs::Vector{S},
                 conv::F=nothing) where {T,S,F}
    if isassigned(p.f)
        p.f[](xs)::S
    else
        execute!(S[], p, xs, conv)
    end
end

retrieve(cs, xs, res, i) =
    isconstant(i) ? cs[constantidx(i)] :
    isinput(i) ? xs[inputidx(i)] :
    res[i.x]

function execute!(res::Vector{S}, p::SLProgram{T}, xs::Vector{S},
                  conv::F=nothing) where {S,T,F}
    # TODO: handle isempty(lines(p))
    empty!(res)

    cs = constants(p)
    for line in lines(p)
        local r::S
        op, i, j = unpack(line)
        x = retrieve(cs, xs, res, i)
        if isexponentiate(op)
            r = x^Int(j.x) # TODO: support bigger j
        elseif isuniplus(op) # serves as assignment (for trivial programs)
            if isa(x, S)
                r = x
            else
                r = conv(x)
            end
        elseif isuniminus(op)
            r = -x
        else
            y = retrieve(cs, xs, res, j)
            if isplus(op)
                r = x + y
            elseif isminus(op)
                r = x - y
            elseif istimes(op)
                r = x * y
            elseif isdivide(op)
                r = divexact(x, y)
            else
                throw(ArgumentError("unknown operation"))
            end
        end
        push!(res, r)
    end
    res[end]
end


## compile!

cretrieve(i) =
    isinput(i) ? Symbol(:x, inputidx(i)) => inputidx(i) :
    isconstant(i) ? Symbol(:c, constantidx(i)) => 0 :
    Symbol(:res, i.x) => 0

# TODO: handle the "conv" argument like in execute!
# (works without, but there can be type-instability)

# return compiled execution function f, and updates
# p.f[] = f, which is not invalidated when p is mutated
function compile!(p::SLProgram)
    res = Expr[]
    fn = :(function (xs::Vector)
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
                :($rk = $x^$(Int(j.x)))
            elseif isuniplus(op)
                :($rk = $x)
            elseif isuniminus(op)
                :($rk = -$x)
            else
                y, idx = cretrieve(j)
                mininput = max(mininput, idx)
                if isplus(op)
                    :($rk = $x + $y)
                elseif isminus(op)
                    :($rk = $x - $y)
                elseif istimes(op)
                    :($rk = $x * $y)
                elseif isdivide(op)
                    :($rk = divexact($x, $y))
                end
            end
        push!(res, line)
    end
    for k = 1:mininput-1
        pushfirst!(res, :($(Symbol(:x, k)) = @inbounds xs[$k]))
    end
    if mininput >= 1
        pushfirst!(res, :($(Symbol(:x, mininput)) = xs[$mininput]))
    end
    append!(fn.args[2].args, res)
    p.f[] = eval(fn)
end
