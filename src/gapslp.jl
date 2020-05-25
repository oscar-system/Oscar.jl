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


## compile!

# compile to an SLProgram

function compile!(gp::GAPSLProgram)
    p = SLProgram()
    local res::Arg

    for i = 1:gp.ngens
        pushop!(p, assign, input(i))
    end

    for line in gp.lines
        if issimpleline(line)
            ptr = p.len + 1
            k = write_list!(p, line)
            res = Arg(ptr)
            if res != k
                pushop!(p, assign, k, res)
            end
            pushop!(p, keep, res)
        elseif isassignline(line)
            list, dst = line
            ptr = p.len
            k = write_list!(p, list)
            res = pushop!(p, assign, k, Arg(dst))
            pushop!(p, keep, Arg(ptr))
        else
            throw(ArgumentError("not implemented"))
        end
    end
    pushfinalize!(p, res)
    gp.slp[] = p
end

function write_list!(p::SLProgram, list::Vector{Int})
    @assert !isempty(list) # TODO: handle empty lists
    n = length(list) >> 1
    local k
    for i = 1:n
        x = Arg(list[2*i-1])
        e = Arg(list[2*i])
        if e.x == 1
            l = x
        else
            l = pushop!(p, exponentiate, x, e)
        end
        k = i == 1 ? l : pushop!(p, times, k, l)
    end
    k
end
