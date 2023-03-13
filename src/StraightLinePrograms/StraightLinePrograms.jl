module StraightLinePrograms

import Base: +, -, *, ^, parent

import ..AbstractAlgebra: evaluate
export SLProgram
export AbstractSLProgram, GAPSLProgram, GAPSLDecision, AtlasSLProgram, AtlasSLDecision, Lazy
export slpcst, slpgen, slpgens, compile, evaluate, nsteps, list, call
import ..@req

abstract type AbstractSLProgram end

(::Type{SLP})(p::AbstractSLProgram) where {SLP<:AbstractSLProgram} =
    compile(SLP, p)

include("straightline.jl")
include("lazy.jl")
include("gap.jl")
include("atlas.jl")


end # module
