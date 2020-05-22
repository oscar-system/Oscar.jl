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
