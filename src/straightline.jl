## constants & predicates

const opmask    = 0xff00000000000000
const argmask   = 0x000000000fffffff
const inputmark = 0x0000000008000000
const tmpmark   = 0x0000000004000000
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
function input(i::Int)
    @assert i < tmpmark
    Arg((i % UInt64) | inputmark)
end

isinput(i::Arg) = inputmark & i.x != 0

# call before mutating, unless p is empty (opposite of pushfinalize!)
function pushinit!(p::SLPoly, k=length(p.cs), koffset=0,
                   idx=eachindex(p.lines))
    offset = first(idx) - 1
    local n, op, i, j
    for outer n in idx
        op, i, j = unpack(p.lines[n])
        if !isinput(i)
            i0 = i.x
            if i0 > k
                i0 -= k-offset
                i0 ⊻= tmpmark
            else
                i0 += koffset
            end
            i = Arg(i0)
        end
        if !isquasiunary(op)
            if !isinput(j)
                j0 = j.x
                if j0 > k
                    j0 -= k-offset
                    j0 ⊻= tmpmark
                else
                    j0 += koffset
                end
                j = Arg(j0)
            end
        end
        p.lines[n] = Line(op, i, j)
    end
    @assert n == lastindex(p.lines)
    if isuniplus(op) && (isinput(i) || (i.x & tmpmark) == 0) # constant
        pop!(p.lines) # discard trivial instruction
        i
    else
        Arg(n % UInt64 | tmpmark)
    end
end

function pushconst!(p::SLPoly{T}, c::T) where T
    push!(p.cs, c)
    l = lastindex(p.cs)
    @assert l < tmpmark
    Arg(l % UInt64)
end

function pushop!(p::SLPoly, op::Op, i::Arg, j::Arg=Arg(0))
    @assert i.x <= argmask && j.x <= argmask
    push!(p.lines, Line(op, i, j))
    l = lastindex(p.lines)
    @assert l < tmpmark
    Arg(l % UInt64 | tmpmark)
end

# make state consistent again
function pushfinalize!(p::SLPoly, ret::Arg)
    k = length(p.cs)
    if isinput(ret) || (ret.x & tmpmark) == 0 ||
            (ret.x ⊻ tmpmark) != lastindex(p.lines)
        # non-trivial return (i.e. no line or not the result of the last line)
        pushop!(p, uniplus, ret)
    end
    for n in eachindex(p.lines)
        op, i, j = unpack(p.lines[n])
        i, j = i.x, j.x
        if i & tmpmark != 0
            i ⊻= tmpmark
            i += k
        end
        if j & tmpmark != 0
            j ⊻= tmpmark
            j += k
        end
        p.lines[n] = pack(op, i, j)
    end
    p
end
