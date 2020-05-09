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
