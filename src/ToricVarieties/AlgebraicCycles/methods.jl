################################################
# 1: Intersections involving closed subvarieties
################################################

function Base.:*(ac::RationalEquivalenceClass, sv::ClosedSubvarietyOfToricVariety)
    if toric_variety(ac) !== toric_variety(sv)
        throw(ArgumentError("The rational equivalence class and the closed subvariety must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ac * RationalEquivalenceClass(sv)
end


function Base.:*(sv::ClosedSubvarietyOfToricVariety, ac::RationalEquivalenceClass)
    if toric_variety(ac) !== toric_variety(sv)
        throw(ArgumentError("The rational equivalence class and the closed subvariety must be defined on the same toric variety, i.e. the same OSCAR variable"))
    end
    return ac * RationalEquivalenceClass(sv)
end
