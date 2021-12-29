#####################
# 1. Defining data of a line bundle
#####################

@doc Markdown.doc"""
    divisor_class(l::ToricLineBundle)

Return the divisor class which defines the toric line bundle `l`.
"""
function divisor_class(l::ToricLineBundle)
    return l.divisor_class
end
export divisor_class

@doc Markdown.doc"""
    variety(l::ToricLineBundle)

Return the toric variety over which the toric line bundle `l` is defined.
"""
function variety(l::ToricLineBundle)
    return l.variety
end
export variety
