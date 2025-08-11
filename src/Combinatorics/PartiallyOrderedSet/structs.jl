abstract type AbstractPartiallyOrderedSet end

mutable struct PartiallyOrderedSet <: AbstractPartiallyOrderedSet
  pm_poset::Polymake.BigObject
  # revidsmap::Vector{Int} # not needed yet
  atomlabels::Vector
  artificial_top::Bool
  artificial_bottom::Bool

  function PartiallyOrderedSet(b::Polymake.BigObject)
    pos = new()
    pos.pm_poset = b
    pos.artificial_top = false
    pos.artificial_bottom = false
    return pos
  end
end

struct PartiallyOrderedSetElement
  parent::PartiallyOrderedSet
  node_id::Int
end
