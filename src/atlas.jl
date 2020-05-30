struct AtlasLine
    cmd::Symbol
    dst::Int
    i::Int
    j::Int
end

const ATLAS_VOID_SLOT = 0

AtlasLine(cmd, dst, i) = AtlasLine(cmd, dst, i, ATLAS_VOID_SLOT)

struct AtlasSLProgram
    code::String # source code
    ngens::Int
    outputs::Vector{Int}
    lines::Vector{AtlasLine}
end

function AtlasSLProgram(code::String)
    codelines = split(code, "\n", keepempty=false)
    labels = String[]
    outputs = Int[]
    lines = AtlasLine[]

    ngens::Int = 0

    function getidx!(label, create=false)
        # no "inp" line, set defaults
        isempty(labels) && push!(labels, "1", "2")
        if ngens == 0
            ngens = length(labels)
        end

        i = findfirst(==(label), labels)
        if i === nothing
            if create
                push!(labels, label)
                lastindex(labels)
            else
                throw(ArgumentError("label $label does not exist"))
            end
        else
            i
        end
    end

    for codeline in codelines
        codeline = split(codeline, '#', limit=2)[1] # remove comments
        codeline = split(codeline, keepempty=false)
        isempty(codeline) && continue
        length(codeline) < 2 && error_invalid_line(codeline)
        cmd = Symbol(codeline[1])

        if cmd == :inp
            ngens != 0 && throw(ArgumentError("\"inp\" line not at the beginning"))
            n = tryparse(Int, codeline[2])
            n === nothing && error_invalid_line(codeline)
            if length(codeline) == 2
                isempty(labels) ||
                    throw(ArgumentError("inp must not omit the names"))
                append!(labels, string.(1:n))
            else
                length(codeline) == 2+n || error_invalid_line(codeline)
                append!(labels, codeline[3:end])
            end
            continue
        end

        if cmd == :oup
            n = tryparse(Int, codeline[2])
            n === nothing && error_invalid_line(codeline)
            if length(codeline) == 2
                isempty(outputs) ||
                    throw(ArgumentError("oup must not omit the names"))
                append!(outputs, getidx!.(string.(1:n)))
            else
                length(codeline) == 2+n || error_invalid_line(codeline)
                append!(outputs, getidx!.(codeline[3:end]))
            end
            continue
        end

        !isempty(outputs) && throw(ArgumentError("\"oup\" line not at the end"))

        line =
            if cmd == :cjr
                check_line_length(codeline, 2)
                args = getidx!.(codeline[2:end], (false, false))
                AtlasLine(:cj, args[1], args[1], args[2])
            elseif cmd in [:cj, :com, :mu]
                check_line_length(codeline, 3)
                args = getidx!.(codeline[2:end], (false, false, true))
                AtlasLine(cmd, args[3], args[1], args[2])
            elseif cmd in [:iv, :cp]
                check_line_length(codeline, 2)
                args = getidx!.(codeline[2:end], (false, true))
                AtlasLine(cmd, args[2], args[1])
            elseif cmd == :pwr
                check_line_length(codeline, 3)
                arg = getidx!(codeline[3])
                dst = getidx!(codeline[4], true)
                AtlasLine(cmd, dst, parse(Int, codeline[2]), arg)
            else
                error_invalid_line(codeline)
            end
        push!(lines, line)
    end

    if isempty(outputs)
        push!(outputs, getidx!("1"), getidx!("2"))
    end
    AtlasSLProgram(code, ngens, outputs, lines)
end

check_line_length(codeline, n) =
    length(codeline) == n+1 || throw(ArgumentError(
        "wrong number of arguments in $(codeline[1]) line"))

error_invalid_line(codeline) =
    throw(ArgumentError("""invalid line: $(join(codeline, " "))"""))


## show

Base.show(io::IO, p::AtlasSLProgram) = print(io, p.code)


## evaluate

evaluate(p::AtlasSLProgram, xs::Vector) = evaluate!(empty(xs), p, xs)

function evaluate!(res::Vector{S}, p::AtlasSLProgram, xs::Vector{S}) where S
    append!(empty!(res), view(xs, 1:p.ngens))

    for l in p.lines
        cmd, dst, i, j = l.cmd, l.dst, l.i, l.j
        a, b = get(res, i, nothing), get(res, j, nothing)

        k =
            if cmd == :cj
                b^-1 * a * b
            elseif cmd == :com
                a^-1 * b^-1 * a * b
            elseif cmd == :iv
                @assert j == ATLAS_VOID_SLOT
                a^-1
            elseif cmd == :mu
                a * b
            elseif cmd == :pwr
                b^i
            elseif cmd == :cp
                @assert j == ATLAS_VOID_SLOT
                a
            else
                @assert false "unexpected command"
            end

        n = lastindex(res)
        if dst == n + 1
            push!(res, k)
        else
            res[dst] = k
        end
    end

    r = 0
    for i in p.outputs
        res[r += 1] = res[i]
    end
    resize!(res, r)
end
