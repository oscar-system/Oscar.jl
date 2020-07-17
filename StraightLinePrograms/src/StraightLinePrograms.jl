module StraightLinePrograms

import Base: +, -, *, ^, parent


using AbstractAlgebra: AbstractAlgebra, Ring, RingElement, RingElem, MPolyRing, MPolyElem, elem_type, Generic

import AbstractAlgebra: base_ring, gen, gens, ngens, nvars, symbols, evaluate, order, addeq!

export LazyPolyRing, LazyPoly, SLPolyRing, SLPoly, SLProgram
export AbstractSLProgram, GAPSLProgram, GAPSLDecision, AtlasSLProgram, AtlasSLDecision, Lazy
export slpcst, slpgen, slpgens, compile, gens, evaluate, nsteps, compose, list, call

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
