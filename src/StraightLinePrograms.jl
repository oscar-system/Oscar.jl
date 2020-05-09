module StraightLinePrograms

import Base: +, -, *, ^, parent

using AbstractAlgebra: Ring, RingElement, MPolyRing, MPolyElem, elem_type

import AbstractAlgebra: base_ring

export LazyPolyRing, LazyPoly


include("lazypolys.jl")

end # module
