abstract type AbstractPartiallyOrderedSet end

mutable struct PartiallyOrderedSet <: AbstractPartiallyOrderedSet
  pm_poset::Polymake.BigObject
  # revidsmap::Vector{Int} # not needed yet
  atomlabels::Vector
  artificial_top::Bool

  function PartiallyOrderedSet(b::Polymake.BigObject)
    pos = new()
    pos.pm_poset = b
    pos.artificial_top = false
    return pos
  end
end

struct PartiallyOrderedSetElement
  parent::PartiallyOrderedSet
  node_id::Int
  # TODO: do we want to store something here?
  #elem::T
end
