module StraightLinePrograms

import Base: +, -, *, ^, parent


using AbstractAlgebra: Ring, RingElement, MPolyRing, MPolyElem, elem_type,
    Generic

import AbstractAlgebra: base_ring, gen, gens, symbols, evaluate, order

export LazyPolyRing, LazyPoly, SLPolyRing, SLPoly, SLProgram
export slpcst, slpgen, slpgens
export GAPSLProgram, GAPSLDecision, AtlasSLProgram, AtlasSLDecision

include("straightline.jl")
include("lazy.jl")
include("lazypolys.jl")
include("slpolys.jl")
include("gap.jl")
include("atlas.jl")


end # module
