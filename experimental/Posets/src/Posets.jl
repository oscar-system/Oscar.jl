module Posets

using ..Oscar

import Base: length, parent
import Oscar: index

export MaximalChainsIterator
export Poset, PosetElem

export poset, poset_elem
export maximal_chains

include("Poset.jl")

end

using .Posets

export MaximalChainsIterator
export Poset, PosetElem

export poset, poset_elem
export maximal_chains
