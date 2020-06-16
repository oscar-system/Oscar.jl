module StraightLinePrograms

import Base: +, -, *, ^, parent


using AbstractAlgebra: Ring, RingElement, MPolyRing, MPolyElem, elem_type,
    Generic

import AbstractAlgebra: base_ring, gen, gens, ngens, nvars, symbols, evaluate, order

export LazyPolyRing, LazyPoly, SLPolyRing, SLPoly, SLProgram
export AbstractSLProgram, GAPSLProgram, GAPSLDecision, AtlasSLProgram, AtlasSLDecision, Free
export slpcst, slpgen, slpgens, compile, gens, evaluate

abstract type AbstractSLProgram end

(::Type{SLP})(p::AbstractSLProgram) where {SLP<:AbstractSLProgram} =
    compile(SLP, p)


include("straightline.jl")
include("lazy.jl")
include("lazypolys.jl")
include("slpolys.jl")
include("gap.jl")
include("atlas.jl")


end # module
