###############################################################################
#
#   LinearAlgebraicGroup, LinearAlgebraicGroupElem
#
###############################################################################

@attributes mutable struct LinearAlgebraicGroup{C<:FieldElem} <: Group
  rs::RootSystem #in the future to be replaced by RootDatum
  G::MatGroup{C} #the actual group
  k::Field #the field over which the group is defined

  # The following fields are not set by default, just for caching
  T::MatGroup{C} #maximal torus
  B::MatGroup{C} #borel subgroup
  U_alphas::Dict{RootSpaceElem,MatGroup} #root subgroups

  function LinearAlgebraicGroup(R::RootSystem, G::MatGroup, k::Field)
    return new{elem_type(k)}(R, G, k)
  end
end

@attributes mutable struct LinearAlgebraicGroupElem{C<:FieldElem} <: GroupElem
  parent::LinearAlgebraicGroup{C}
  mat::MatGroupElem{C} #the actual element

  function LinearAlgebraicGroupElem(
    parent::LinearAlgebraicGroup{C}, MGE::MatGroupElem{C}; check::Bool=true
  ) where {C<:FieldElem}
    if check
      @req MGE in parent.G "The given matrix group element is not an element of the group"
    end
    return new{C}(parent, MGE)
  end
end
