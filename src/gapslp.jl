const ReturnList = Vector{Vector{Int}}

const GAPStraightLine = Union{Vector{Int},            # e.g. [1, 2, 2, -1]
                              Tuple{Vector{Int},Int}, # e.g. ([1, 2], 2)
                              ReturnList}             # return list

issimpleline(line) = line isa Vector{Int}
isassignline(line) = line isa Tuple
isreturnline(line) = line isa ReturnList

struct GAPSLProgram
    lines::Vector{GAPStraightLine}
    ngens::Int
    slp::Ref{SLProgram}
end

function GAPSLProgram(lines::Vector, ngens::Integer=-1)
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
        if ngens < 0 && i < n && l isa Vector{Int}
            throw(ArgumentError("ngens must be specified"))
        end
        if l isa Tuple
            ng = max(ng, maxnohave(l[1]))
            if ngens < 0
                union!(have, 1:ng)
            elseif ng > 0
                undef_slot()
            end
            push!(have, l[2])
        elseif l isa Vector{Int}
            ng = max(ng, maxnohave(l))
            ngens >= 0 && ng > 0 && undef_slot()
            push!(have, maximum(have)+1)
        else # ReturnList
            for li in l
                ng = max(ng, maxnohave(li))
                ngens >= 0 && ng > 0 && undef_slot()
            end
        end
    end
    GAPSLProgram(ls, ngens < 0 ? ng : ngens, Ref{SLProgram}())
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

function pushline!(lines::Vector{GAPStraightLine}, line::Vector)
    length(line) == 2 && line[1] isa Vector{Int} && line[2] isa Int ||
        invalid_list_error(line)
    pushline!(lines, (line[1], line[2]))
end


## show

expshow(x, i) = i == 1 ? "r[$x]" : "r[$x]^$i"
prodshow(io, x) = join(io, (expshow(x[2*i-1], x[2*i]) for i=1:(length(x)>>1)), '*')

function Base.show(io::IO, p::GAPSLProgram)
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
            l, dst = line
            # TODO: when l > k
            print(io, "\nr[$dst] = ")
            prodshow(io, l)
        else
            print("\n# return values:\n[ ")
            for (i, l) in enumerate(line)
                i == 1 || print(io, ", ")
                prodshow(io, l)
            end
            return print(io, " ]")
        end
    end
    print("\n# return value:\nr[$k]")
end


## evaluate

evaluate(p::GAPSLProgram, xs::Vector) = evaluate!(empty(xs), p, xs)
# TODO: use compiled version if present?

expterm(x, i) = i == 1 ? x : x^i
prodlist(res, x) = prod(expterm(res[x[2*i-1]], x[2*i]) for i=1:(length(x)>>1))

function evaluate!(res::Vector{S}, p::GAPSLProgram, xs::Vector{S}) where {S}
    empty!(res)

    for i = 1:p.ngens
        push!(res, xs[i])
    end

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
    return res[k]
end


## compile!

# compile to an SLProgram

function compile!(gp::GAPSLProgram)
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
    multi || pushfinalize!(p, res)
    gp.slp[] = p
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
