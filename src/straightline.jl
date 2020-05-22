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


## building SLPoly

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

function pushinit!(p::SLPoly)
    op, i, j = unpack(p.lines[end])
    if isuniplus(op) && (isinput(i) || isconstant(i))
        pop!(p.lines) # discard trivial instruction
        i
    else
        Arg(lastindex(p.lines) % UInt64)
    end
end

function combine!(p::SLPoly{T}, q::SLPoly{T}) where T
    i1 = pushinit!(p)
    koffset = length(p.cs)
    len = length(p.lines)
    append!(p.lines, q.lines)
    append!(p.cs, q.cs)

    @assert length(p.cs) < cstmark # TODO: should not be @assert
    for n = len+1:lastindex(p.lines)
        op, i, j = unpack(p.lines[n])
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
        p.lines[n] = Line(op, i, j)
        # TODO: write conditionally only when modifications
    end
    i2 = pushinit!(p)
    i1, i2
end

function pushconst!(p::SLPoly{T}, c::T) where T
    push!(p.cs, c)
    l = lastindex(p.cs)
    @assert l < cstmark
    asconstant(l)
end

function pushop!(p::SLPoly, op::Op, i::Arg, j::Arg=Arg(0))
    @assert i.x <= argmask && j.x <= argmask
    push!(p.lines, Line(op, i, j))
    l = lastindex(p.lines)
    @assert l < cstmark
    Arg(l % UInt64)
end

# make state consistent again
function pushfinalize!(p::SLPoly, ret::Arg)
    k = length(p.cs)
    if isinput(ret) || isconstant(ret) || ret.x != lastindex(p.lines)
        # non-trivial return (i.e. no line or not the result of the last line)
        pushop!(p, uniplus, ret)
    end
    p
end
