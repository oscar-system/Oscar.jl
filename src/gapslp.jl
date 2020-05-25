const ReturnList = Vector{Vector{Int}}

const GAPStraightLine = Union{Vector{Int},            # e.g. [1, 2, 2, -1]
                              Tuple{Vector{Int},Int}, # e.g. ([1, 2], 2)
                              ReturnList}             # return list

struct GAPSLProgram
    lines::Vector{GAPStraightLine}
    ngens::Int
    slp::Ref{SLProgram}
end

function GAPSLProgram(lines::Vector, ngens::Integer=-1)
    ls = GAPStraightLine[]
    n = length(lines)
    ng = 0
    if ngens < 0
        have = BitSet()
    end

    function maxnohave(word)
        ng = 0
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

    for (i, line) in enumerate(lines)
        pushline!(ls, line)
        ngens < 0 || continue
        l = ls[end]
        if i < n && l isa Vector{Int}
            throw(ArgumentError("ngens must be specified"))
        end
        if l isa Tuple
            ng = max(ng, maxnohave(l[1]))
            union!(have, 1:ng)
            push!(have, l[2])
        elseif l isa Vector{Int}
            ng = max(ng, maxnohave(l))
        else # ReturnList
            for li in l
                ng = max(ng, maxnohave(li))
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
