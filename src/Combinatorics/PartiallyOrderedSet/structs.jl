abstract type AbstractPartiallyOrderedSet end

struct PartiallyOrderedSet <: AbstractPartiallyOrderedSet
  pm_poset::Polymake.BigObject
end

#=
do we need this?
struct PartiallyOrderedSetElement{T}
  parent::PartiallyOrderedSet
  node_id::Int
  elem::T
end
=#
