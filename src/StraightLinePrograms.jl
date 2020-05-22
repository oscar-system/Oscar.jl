module StraightLinePrograms

import Base: +, -, *, ^, parent

using AbstractAlgebra: Ring, RingElement, MPolyRing, MPolyElem, elem_type,
    Generic

import AbstractAlgebra: base_ring, gen, gens, symbols, evaluate

export LazyPolyRing, LazyPoly, SLPolyRing, SLPoly


include("recpolys.jl")
include("lazypolys.jl")
include("slprograms.jl")
include("slpolys.jl")
include("straightline.jl")

end # module
