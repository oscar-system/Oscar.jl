module StraightLinePrograms

import Base: +, -, *, ^, parent

import ..AbstractAlgebra: evaluate, add!, sub!, mul!, neg!

export AbstractSLProgram
export AtlasSLDecision
export AtlasSLProgram
export GAPSLDecision
export GAPSLProgram
export Lazy
export SLProgram
export call
export compile
export evaluate
export list
export nsteps
export slpcst
export slpgen
export slpgens

import ..@req

abstract type AbstractSLProgram end

(::Type{SLP})(p::AbstractSLProgram) where {SLP<:AbstractSLProgram} =
    compile(SLP, p)

include("straightline.jl")
include("lazy.jl")
include("gap.jl")
include("atlas.jl")


end # module
