################################################
# 1: Intersections involving closed subvarieties
################################################

function Base.:*(ac::RationalEquivalenceClass, sv::ClosedSubvarietyOfToricVariety)
  @req toric_variety(ac) === toric_variety(sv) "The rational equivalence class and the closed subvariety must be defined on the same toric variety"
  return ac * rational_equivalence_class(sv)
end

function Base.:*(sv::ClosedSubvarietyOfToricVariety, ac::RationalEquivalenceClass)
  @req toric_variety(ac) === toric_variety(sv) "The rational equivalence class and the closed subvariety must be defined on the same toric variety"
  return ac * rational_equivalence_class(sv)
end
