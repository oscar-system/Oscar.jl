module StraightLinePrograms

import Base: +, -, *, ^, parent

using AbstractAlgebra: Ring, RingElement, MPolyRing, MPolyElem, elem_type

import AbstractAlgebra: base_ring, gen, gens, symbols

export LazyPolyRing, LazyPoly, SLPolyRing, SLPoly


include("lazypolys.jl")
include("slpolys.jl")
include("straightline.jl")

end # module
