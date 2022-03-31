module StraightLinePrograms

import Base: +, -, *, ^, parent

export SLProgram
export AbstractSLProgram, GAPSLProgram, GAPSLDecision, AtlasSLProgram, AtlasSLDecision, Lazy
export slpcst, slpgen, slpgens, compile, evaluate, nsteps, list, call

abstract type AbstractSLProgram end

(::Type{SLP})(p::AbstractSLProgram) where {SLP<:AbstractSLProgram} =
    compile(SLP, p)

include("straightline.jl")
include("lazy.jl")
include("gap.jl")
include("atlas.jl")


end # module
