@doc Markdown.doc"""
    adding(x::Int,y::Int)

Return the sum of two integers.

```jldoctest
julia> adding(1,2)
3
```
"""
function adding(x::Int,y::Int)
    return x + y
end
export adding
