## SLPolyRing (SL = straight-line)

struct SLPolyRing{T<:RingElement,R<:Ring} <: MPolyRing{T}
    base_ring::R
    S::Vector{Symbol}

    SLPolyRing(r::Ring, s::Vector{Symbol}) = new{elem_type(r),typeof(r)}(r, s)
end

base_ring(S::SLPolyRing) = S.base_ring

symbols(S::SLPolyRing) = S.S

(PR::SLPolyRing{T})() where {T} = SLPoly(PR, T[], UInt64[])


## SLPoly

struct SLPoly{T<:RingElement,SLPR<:SLPolyRing{T}} <: MPolyElem{T}
    parent::SLPR
    cs::Vector{T}          # constants
    lines::Vector{UInt64}  # instructions
    f::Ref{Function}       # compiled evalutation

    SLPoly(parent, xs, c) =
        new{elem_type(base_ring(parent)),typeof(parent)}(
            parent, xs, c, Ref{Function}())
end

parent(p::SLPoly) = p.parent

function check_parent(p::SLPoly, q::SLPoly)
    p.parent === q.parent ||
        throw(ArgumentError("incompatible parents"))
    p.parent
end

function Base.copy!(p::SLPoly{T}, q::SLPoly{T}) where {T}
    check_parent(p, q)
    copy!(p.cs, q.cs)
    copy!(p.lines, q.lines)
    p
end

function Base.copy(q::SLPoly)
    p = q.parent()
    copy!(p, q)
    p
end


## show

function Base.show(io::IO, ::MIME"text/plain", p::SLPoly{T}) where T
    n = length(p.lines)
    R = parent(p)
    syms = symbols(R)

    # pre-compute line representations via LazyPoly
    L = LazyPolyRing(base_ring(R))
    xs = map(L, syms)
    res = empty(xs)
    evaluate!(res, p, xs, L)

    for (k, line) in enumerate(p.lines)
        sk = string(k)
        print(io, ' '^max(0, 3-length(sk)), '#', sk, " = ")
        op, i, j = unpack(line)
        x = showarg(p.cs, syms, i)
        y = isunary(op) ? "" :
            isquasiunary(op) ? string(j) :
            showarg(p.cs, syms, j)
        print(io, showop[op], ' ', x, ' ', y)
        print(io, "\t==>\t", res[k+length(p.cs)])
        k == n || println(io)
    end
end

function showarg(cs, syms, i)
    n = length(cs)
    if i <= n
        string(cs[i])
    elseif i & inputmark == 0
        "#$(i-n)"
    else
        string(syms[i ⊻ inputmark])
    end
end


## mutating ops

function combine!(p::SLPoly{T}, q::SLPoly{T}) where T
    i = pushinit!(p)
    koffset = length(p.cs)
    len = length(p.lines)
    append!(p.lines, q.lines)
    append!(p.cs, q.cs)
    j = pushinit!(p, length(q.cs), koffset, len+1:lastindex(p.lines))
    i, j
end

function addeq!(p::SLPoly{T}, q::SLPoly{T}) where {T}
    pushop!(p, plus, combine!(p, q)...)
    pushfinalize!(p)
end

function subeq!(p::SLPoly{T}, q::SLPoly{T}) where {T}
    pushop!(p, minus, combine!(p, q)...)
    pushfinalize!(p)
end

function subeq!(p::SLPoly)
    i = pushinit!(p)
    pushop!(p, uniminus, i)
    pushfinalize!(p)
end

function muleq!(p::SLPoly{T}, q::SLPoly{T}) where {T}
    pushop!(p, times, combine!(p, q)...)
    pushfinalize!(p)
end

function expeq!(p::SLPoly, e::Integer)
    i = pushinit!(p)
    pushop!(p, exponentiate, i, Int(e))
    pushfinalize!(p)
end


## unary/binary ops

+(p::SLPoly{T}, q::SLPoly{T}) where {T} = addeq!(copy(p), q)

*(p::SLPoly{T}, q::SLPoly{T}) where {T} = muleq!(copy(p), q)

-(p::SLPoly{T}, q::SLPoly{T}) where {T} = subeq!(copy(p), q)

-(p::SLPoly) = subeq!(copy(p))

^(p::SLPoly, e::Integer) = expeq!(copy(p), e)


## evaluate

retrieve(xs, res, i) =
    (i & inputmark) == 0 ? res[i] : xs[i ⊻ inputmark]

function evaluate!(res::Vector{S}, p::SLPoly{T}, xs::Vector{S},
                   R::Ring=parent(xs[1])) where {S,T}
    # TODO: handle isempty(p.lines)
    resize!(res, length(p.cs))
    for i in eachindex(res)
        res[i] = R(p.cs[i])
    end

    for line in p.lines
        op, i, j = unpack(line)
        x = retrieve(xs, res, i)
        if isexponentiate(op)
            r = x^Int(j) # TODO: support bigger j
        elseif isuniplus(op) # serves as assignment (for trivial programs)
            r = x
        elseif isuniminus(op)
            r = -x
        else
            y = retrieve(xs, res, j)
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


## conversion Lazy -> SLP

function (R::SLPolyRing{T})(p::LazyPoly{T}) where {T}
    q = R()
    pushpoly!(q, p.p)
    pushfinalize!(q)
end

pushpoly!(q, p::Const) = pushconst!(q, p.c)

pushpoly!(q, p::Gen) = input(findfirst(==(p.g), symbols(parent(q))))

function pushpoly!(q, p::Union{PlusPoly,TimesPoly})
    # TODO: handle isempty(p.xs) ?
    op = p isa PlusPoly ? plus : times
    x, xs = Iterators.peel(p.xs)
    i = pushpoly!(q, x)
    for x in xs
        j = pushpoly!(q, x)
        i = pushop!(q, op, i, j)
    end
    i
end

function pushpoly!(q, p::MinusPoly)
    i = pushpoly!(q, p.p)
    j = pushpoly!(q, p.q)
    pushop!(q, minus, i, j)
end

pushpoly!(q, p::UniMinusPoly) = pushop!(q, uniminus, pushpoly!(q, p.p))

pushpoly!(q, p::ExpPoly) = pushop!(q, exponentiate, pushpoly!(q, p.p), p.e)


## conversion SLPoly -> MPoly

function Base.convert(R::MPolyRing, p::SLPoly)
    symbols(R) == symbols(parent(p)) ||
        throw(ArgumentError("incompatible symbols"))
    xs = gens(R)
    res = empty(xs)
    evaluate!(res, p, xs, R)
end
