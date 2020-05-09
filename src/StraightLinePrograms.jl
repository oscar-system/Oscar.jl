module StraightLinePrograms

import Base: +, -, *, ^

using AbstractAlgebra: Ring, RingElement, MPolyRing, elem_type

import AbstractAlgebra: base_ring

export LazyPolyRing


include("lazypolys.jl")

end # module
