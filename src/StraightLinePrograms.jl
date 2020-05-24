module StraightLinePrograms

import Base: +, -, *, ^, parent


using AbstractAlgebra: Ring, RingElement, MPolyRing, MPolyElem, elem_type,
    Generic

import AbstractAlgebra: base_ring, gen, gens, symbols, evaluate

export LazyPolyRing, LazyPoly, SLPolyRing, SLPoly, SLProgram
export slpcst, slpgen, slpgens


include("lazy.jl")
include("straightline.jl")
include("lazypolys.jl")
include("slpolys.jl")


end # module
