## Op, Line, Arg

struct Op
    x::UInt64
end

struct Line
    x::UInt64
end

struct Arg
    x::UInt64
end


struct SLProgram{T}
    cs::Vector{T}       # constants
    lines::Vector{Line} # instructions
    f::Ref{Function}    # compiled execution
end

SLProgram(cs, lines) = SLProgram(cs, lines, Ref{Function}())

constants(p::SLProgram) = p.cs
lines(p::SLProgram) = p.lines
