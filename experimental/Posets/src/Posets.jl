module Posets

using ..Oscar

import Base: length, parent
import Oscar: index

export MaximalChainsIterator
export Poset, PosetElem

export maximal_chains
export poset

include("Poset.jl")

end

using .Posets

export MaximalChainsIterator
export Poset, PosetElem

export maximal_chains
export poset
