const ReturnList = Vector{Vector{Int}}

const GAPStraightLine = Union{Vector{Int},            # e.g. [1, 2, 2, -1]
                              Tuple{Vector{Int},Int}, # e.g. ([1, 2], 2)
                              Tuple{Int,Int},         # (i, n) for ["Order", i, n]
                              ReturnList}             # return list

issimpleline(line) = line isa Vector{Int}
isassignline(line) = line isa Tuple{Vector{Int},Int}
isreturnline(line) = line isa ReturnList
isorderline(line)  = line isa Tuple{Int,Int}

abstract type AbstractGAPSL <: AbstractSLProgram end

struct GAPSLProgram <: AbstractGAPSL
    lines::Vector{GAPStraightLine}
    ngens::Int
    slp::Ref{SLProgram}
end

struct GAPSLDecision <: AbstractGAPSL
    lines::Vector{GAPStraightLine}
    ngens::Int
    slp::Ref{SLProgram}
end

function (::Type{SLP})(lines::Vector, ngens::Integer=-1) where {SLP<:AbstractGAPSL}
    ls = GAPStraightLine[]
    n = length(lines)
    ng = 0

    have = ngens < 0 ? BitSet(0) : BitSet(0:ngens)

    function maxnohave(word)
        local ng = 0
        skip = true
        for i in word
            skip = !skip
            skip && continue
            if !(i in have)
                ng = max(ng, i)
            end
        end
        ng
    end

    undef_slot() = throw(ArgumentError("invalid use of undef memory"))

    for (i, line) in enumerate(lines)
        pushline!(ls, line)
        l = ls[end]
        if ngens < 0 && issimpleline(l) && (i < n || SLP === GAPSLDecision)
            throw(ArgumentError("ngens must be specified"))
        end
        if isassignline(l)
            ng = max(ng, maxnohave(l[1]))
            if ngens < 0
                union!(have, 1:ng)
            elseif ng > 0
                undef_slot()
            end
            push!(have, l[2])
        elseif issimpleline(l)
            ng = max(ng, maxnohave(l))
            ngens >= 0 && ng > 0 && undef_slot()
            push!(have, maximum(have)+1)
        elseif isorderline(l)
            SLP <: GAPSLDecision ||
                throw(ArgumentError("\"Order\" line only allowed in GAPSLDecision"))
            # GAP doesn't seem to take an "Order" line into account for determining
            # the number of generators
        else
            @assert isreturnline(l) "unknown line"
            SLP <: GAPSLProgram ||
                throw(ArgumentError("return list only allowed in GAPSLProgram"))
            for li in l
                ng = max(ng, maxnohave(li))
                ngens >= 0 && ng > 0 && undef_slot()
            end
        end
    end
    SLP(ls, ngens < 0 ? ng : ngens, Ref{SLProgram}())
end

invalid_list_error(line) = throw(ArgumentError("invalid line or list: $line"))

check_element(list) = iseven(length(list)) || invalid_list_error(list)

check_last(lines) = isempty(lines) || !isa(lines[end], ReturnList) ||
    throw(ArgumentError("return list only allowed in last position"))

function pushline!(lines::Vector{GAPStraightLine}, line::Vector{Int})
    check_last(lines)
    check_element(line)
    push!(lines, line)
end

function pushline!(lines::Vector{GAPStraightLine}, line::Vector{Vector{Int}})
    check_last(lines)
    for l in line
        check_element(l)
    end
    push!(lines, line)
end

function pushline!(lines::Vector{GAPStraightLine}, line::Tuple{Vector{Int},Int})
    check_last(lines)
    check_element(line[1])
    push!(lines, line)
end

function pushline!(lines::Vector{GAPStraightLine}, line::Tuple{Int,Int})
    check_last(lines) # redundant
    push!(lines, line)
end

# catch-all function to parse input in GAP style
function pushline!(lines::Vector{GAPStraightLine}, line::AbstractVector)
    newline =
        if length(line) == 2 &&
                line[1] isa AbstractVector && line[2] isa Integer
            all(isinteger, line[1]) || invalid_list_error(line)
            (Vector{Int}(line[1]), Int(line[2])) # assign line
        elseif length(line) == 3 && line[1] == "Order"
            isinteger(line[2]) && isinteger(line[3]) ||
                invalid_list_error(line)
            (Int(line[2]), Int(line[3]))
        elseif isempty(line) || line[1] isa AbstractVector # return line
            all(v -> v isa AbstractVector && all(isinteger, v), line) ||
                invalid_list_error(line)
            Vector{Int}[Vector{Int}(l) for l in line]
        elseif line isa AbstractVector # simple line
            all(isinteger, line) || invalid_list_error(line)
            Vector{Int}(line)
        else
            invalid_list_error(line)
        end
    pushline!(lines, newline)
end


## show

expshow(x, i) = i == 1 ? "r[$x]" : "r[$x]^$i"
prodshow(io, x) = join(io, (expshow(x[2*i-1], x[2*i]) for i=1:(length(x)>>1)), '*')

function Base.show(io::IO, p::AbstractGAPSL)
    k = p.ngens
    println(io, "# input:")
    print(io, "r = [ ")
    join(io, ("g$i" for i in 1:k), ", ")
    println(io, " ]")
    print("# program:")

    for line in p.lines
        if issimpleline(line)
            k += 1
            print(io, "\nr[$k] = ")
            prodshow(io, line)
        elseif isassignline(line)
            l, k = line
            print(io, "\nr[$k] = ")
            prodshow(io, l)
        elseif isorderline(line)
            x, e = line
            print(io, "\norder( r[$x] ) == $e || return false")
        elseif isreturnline(line)
            print("\n# return values:\n[ ")
            for (i, l) in enumerate(line)
                i == 1 || print(io, ", ")
                prodshow(io, l)
            end
            return print(io, " ]")
        else
            @assert false "unknown line"
        end
    end
    if p isa GAPSLProgram
        print("\n# return value:\nr[$k]")
    else
        print("\n# return value:\ntrue")
    end
end


## evaluate

evaluate(p::AbstractGAPSL, xs::Vector) = evaluate!(empty(xs), p, xs)
# TODO: use compiled version if present?

expterm(x, i) = i == 1 ? x : x^i
prodlist(res, x) = prod(expterm(res[x[2*i-1]], x[2*i]) for i=1:(length(x)>>1))

function evaluate!(res::Vector{S}, p::AbstractGAPSL, xs::Vector{S}) where {S}
    append!(empty!(res), view(xs, 1:p.ngens))

    local k
    for line in p.lines
        if issimpleline(line)
            push!(res, prodlist(res, line))
            k = lastindex(res)
        elseif isassignline(line)
            k = line[2]
            r = prodlist(res, line[1])
            if k == lastindex(res) + 1
                push!(res, r)
            else
                res[k] = r
            end
        elseif isorderline(line)
            x = res[line[1]]
            order(x) == line[2] || return false
        else
            k = lastindex(res) + 1
            for l in line
                push!(res, prodlist(res, l))
            end
            if k > 1
                for i in 1:length(line)
                    res[i] = res[k]
                    k += 1
                end
            end
            resize!(res, length(line))
            return res
        end
    end
    if p isa GAPSLProgram
        res[k]
    else
        true
    end
end


## compile

# compile to an SLProgram

compile!(gp::AbstractGAPSL) = gp.slp[] = compile(gp)

compile(gp::AbstractGAPSL) = compile(SLProgram, gp)

function compile(::Type{SLProgram}, gp::AbstractGAPSL)
    p = SLProgram()
    local res::Arg

    for i = 1:gp.ngens
        pushop!(p, assign, input(i))
    end

    multi = false
    for line in gp.lines
        @assert !multi
        if issimpleline(line)
            res = Arg(p.len + 1)
            k = write_list!(p, line)
            if res != k
                pushop!(p, assign, k, res)
            end
            pushop!(p, keep, res)
        elseif isassignline(line)
            list, dst = line
            res = Arg(dst)
            ptr = Arg(max(p.len, dst))
            k = write_list!(p, list)
            if res != k
                pushop!(p, assign, k, res)
            end
            pushop!(p, keep, ptr)
        elseif isorderline(line)
            x = Arg(line[1])
            ord = pushint!(p, line[2])
            pushop!(p, decision, x, ord)
        else
            reslist = Arg[]
            for l in line
                k = write_list!(p, l)
                push!(reslist, k)
            end
            for (i, r) in enumerate(reslist)
                @assert i != r.x # can happen? if so, don't assign
                pushop!(p, assign, r, Arg(i))
            end
            pushop!(p, keep, Arg(length(reslist)))
            multi = true
        end
    end
    if gp isa GAPSLDecision
        setdecision!(p)
    elseif !multi
        pushfinalize!(p, res)
    end
    p
end

function write_list!(p::SLProgram, list::Vector{Int})
    @assert !isempty(list) # TODO: handle empty lists
    n = length(list) >> 1
    local k
    for i = 1:n
        x = Arg(list[2*i-1])
        e = list[2*i]
        if e == 1
            l = x
        else
            l = pushop!(p, exponentiate, x, intarg(e))
        end
        k = i == 1 ? l : pushop!(p, times, k, l)
    end
    k
end
