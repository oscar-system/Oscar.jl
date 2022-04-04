@doc Markdown.doc"""
    contains(tvs::ToricVanishingSet, l::ToricLineBundle)

Checks if the toric line bundle `l` is contained in the toric vanishing set `tvs`.
"""
function contains(tvs::ToricVanishingSet, l::ToricLineBundle)
    class = divisor_class(l).coeff
    class = [class[1,i] for i in 1:ncols(class)]
    for p in polyhedra(tvs)
        if contains(p, class)
            return false
        end
    end
    return true
end
export contains