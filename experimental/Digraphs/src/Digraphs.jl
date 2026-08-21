using Oscar: AbstractGraph, Directed, GAPGroup, GAPGroupElem,
             Graph, MixedGraph, Mixed, Polymake, Undirected, pm_object, @req

import GAP
import GAP: GapObj

include("DigraphWrap.jl")
include("Types.jl")
include("Constructors.jl")
include("Operators.jl")
include("Attributes.jl")
include("Properties.jl")
include("Isomorphisms.jl")
include("Export.jl")
