## SLPolyRing (SL = straight-line)

struct SLPolyRing{T<:RingElement,R<:Ring} <: MPolyRing{T}
    base_ring::R
    S::Vector{Symbol}

    SLPolyRing(r::Ring, s::Vector{Symbol}) = new{elem_type(r),typeof(r)}(r, s)
end

base_ring(S::SLPolyRing) = S.base_ring

symbols(S::SLPolyRing) = S.S

(S::SLPolyRing{T})(c::T=zero(base_ring(S))) where {T} = S(Const(c))

function gen(S::SLPolyRing{T}, i::Integer) where {T}
    s = symbols(S)[i]
    S(Gen{T}(s))
end

gens(S::SLPolyRing{T}) where {T} = [S(Gen{T}(s)) for s in symbols(S)]


## SLPoly

struct SLPoly{T<:RingElement,SLPR<:SLPolyRing{T}} <: MPolyElem{T}
    parent::SLPR
    cs::Vector{T}          # constants
    lines::Vector{UInt64}  # instructions
    f::Ref{Function}       # compiled evalutation

    SLPoly(parent, cs, lines) =
        new{elem_type(base_ring(parent)),typeof(parent)}(
            parent, cs, lines, Ref{Function}())
end

# create invalid poly
SLPoly(parent::SLPolyRing{T}) where {T} = SLPoly(parent, T[], UInt64[])

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
    p = SLPoly(q.parent)
    copy!(p, q)
    p
end


## show

function evaluate_lazy(p)
    # pre-compute line representations via LazyPoly
    R = parent(p)
    L = LazyPolyRing(base_ring(R))
    xs = map(L, symbols(R))
    res = empty(xs)
    evaluate!(res, p, xs, L)
    res
end

Base.show(io::IO, p::SLPoly) = show(io, evaluate_lazy(p)[end])

function Base.show(io::IO, ::MIME"text/plain", p::SLPoly{T}) where T
    n = length(p.lines)
    syms = symbols(parent(p))
    res = evaluate_lazy(p)
    if n == 1
        # trivial program, show only result
        return show(io, res[end])
    end

    for (k, line) in enumerate(p.lines)
        op, i, j = unpack(line)
        if 1 < k == n && op == uniplus && i == k-1+length(p.cs)
            # 1 < k for trivial SLPs returning a constant
            break
        end
        k == 1 || println(io)
        sk = string(k)
        print(io, ' '^max(0, 3-length(sk)), '#', sk, " = ")
        x = showarg(p.cs, syms, i)
        y = isunary(op) ? "" :
            isquasiunary(op) ? string(j) :
            showarg(p.cs, syms, j)
        print(io, showop[op], ' ', x, ' ', y)
        print(io, "\t==>\t", res[k+length(p.cs)])
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
    i = pushop!(p, plus, combine!(p, q)...)
    pushfinalize!(p, i)
end

function subeq!(p::SLPoly{T}, q::SLPoly{T}) where {T}
    i = pushop!(p, minus, combine!(p, q)...)
    pushfinalize!(p, i)
end

function subeq!(p::SLPoly)
    i = pushinit!(p)
    i = pushop!(p, uniminus, i)
    pushfinalize!(p, i)
end

function muleq!(p::SLPoly{T}, q::SLPoly{T}) where {T}
    i = pushop!(p, times, combine!(p, q)...)
    pushfinalize!(p, i)
end

function expeq!(p::SLPoly, e::Integer)
    i = pushinit!(p)
    i = pushop!(p, exponentiate, i, Int(e))
    pushfinalize!(p, i)
end


## unary/binary ops

+(p::SLPoly{T}, q::SLPoly{T}) where {T} = addeq!(copy(p), q)

*(p::SLPoly{T}, q::SLPoly{T}) where {T} = muleq!(copy(p), q)

-(p::SLPoly{T}, q::SLPoly{T}) where {T} = subeq!(copy(p), q)

-(p::SLPoly) = subeq!(copy(p))

^(p::SLPoly, e::Integer) = expeq!(copy(p), e)


## adhoc ops

+(p::SLPoly{T}, q::T) where {T<:RingElement} = p + parent(p)(q)
+(q::T, p::SLPoly{T}) where {T<:RingElement} = parent(p)(q) + p

-(p::SLPoly{T}, q::T) where {T<:RingElement} = p - parent(p)(q)
-(q::T, p::SLPoly{T}) where {T<:RingElement} = parent(p)(q) - p

*(p::SLPoly{T}, q::T) where {T<:RingElement} = p * parent(p)(q)
*(q::T, p::SLPoly{T}) where {T<:RingElement} = parent(p)(q) * p


## evaluate

function evaluate(p::SLPoly{T}, xs::Vector{T}) where {T <: RingElement}
    if isassigned(p.f)
        p.f[](xs)::T
    else
        evaluate!(T[], p, xs)
    end
end

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

(R::SLPolyRing{T})(p::LazyPoly{T}) where {T} = R(p.p)

function (R::SLPolyRing{T})(p::RecPoly{T}) where {T}
    q = SLPoly(R)
    i = pushpoly!(q, p)
    pushfinalize!(q, i)
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


## conversion MPoly -> SLPoly

function Base.convert(R::SLPolyRing, p::Generic.MPoly)
    # TODO: currently handles only default ordering
    symbols(R) == symbols(parent(p)) ||
        throw(ArgumentError("incompatible symbols"))
    q = SLPoly(R)
    @assert lastindex(p.coeffs) < tmpmark
    @assert isempty(q.cs)
    copy!(q.cs, p.coeffs)
    @assert length(q.cs) == size(p.coeffs, 1)
    exps = UInt64[]
    monoms = [Pair{UInt64,UInt64}[] for _ in axes(p.exps, 1)]
    for v in reverse(axes(p.exps, 1))
        copy!(exps, view(p.exps, v, :))
        unique!(sort!(exps))
        for e in exps
            e == 0 && continue
            k = pushop!(q, exponentiate, input(size(p.exps, 1) + 1 - v), e)
            push!(monoms[v], e => k)
        end
        # TODO: handle constants
    end
    k = 0
    for t in eachindex(q.cs)
        i = t
        j = 0
        for v in reverse(axes(p.exps, 1))
            e = p.exps[v, t]
            if  e != 0
                j = monoms[v][searchsortedfirst(monoms[v], e, by=first)][2]
                i = pushop!(q, times, i, j)
            end
        end
        if k == 0
            k = i
        else
            k = pushop!(q, plus, k, i)
        end
    end
    pushfinalize!(q, k)
    q
end


## conversion SLPoly -> MPoly

function Base.convert(R::MPolyRing, p::SLPoly)
    symbols(R) == symbols(parent(p)) ||
        throw(ArgumentError("incompatible symbols"))
    xs = gens(R)
    res = empty(xs)
    evaluate!(res, p, xs, R)
end


## compile!

cretrieve(i) = i & inputmark == 0 ?
    Symbol(:res, i) => 0 :
    Symbol(:x, i ⊻ inputmark) => i ⊻ inputmark

# return compiled evaluation function f, and updates
# p.f[] = f, which is not invalidated when p is mutated
function compile!(p::SLPoly)
    res = Expr[]
    fn = :(function (xs::Vector)
           end)
    k = 0
    cs = p.cs
    for c in p.cs
        k += 1
        push!(res, :($(Symbol(:res, k)) = $cs[$k]))
    end
    mininput = 0
    for line in p.lines
        k += 1
        rk = Symbol(:res, k)
        op, i, j = unpack(line)
        x, idx = cretrieve(i)
        mininput = max(mininput, idx)
        line =
            if isexponentiate(op)
                :($rk = $x^$(Int(j)))
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
