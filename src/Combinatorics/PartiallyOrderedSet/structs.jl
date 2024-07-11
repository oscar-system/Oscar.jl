abstract type AbstractPartiallyOrderedSet end

struct PartiallyOrderedSet <: AbstractPartiallyOrderedSet
  pm_poset::Polymake.BigObject
end

struct PartiallyOrderedSetElement{T}
  parent::PartiallyOrderedSet
  node_id::Int
  elem::T
end
