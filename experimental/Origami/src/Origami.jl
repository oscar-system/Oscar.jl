# Add your new types, functions, and methods here.

export origami, veech_group

module OrigamiHelper

using ..GAP

function __init__()
  GAP.Packages.install("https://ag-weitze-schmithusen.github.io/ModularGroup/PackageInfo.g")
  GAP.Packages.install("https://ag-weitze-schmithusen.github.io/Origami/PackageInfo.g")
  GAP.Packages.load("Origami")
end

end

include("types.jl")

@doc raw"""
    origami(h::PermGroupElem, v::PermGroupElem)

blablabla
# Examples
```jldoctest
julia> h = @perm (1,2)
(1,2)

julia> v = @perm (1,2)
(1,2)

julia> o = origami(h,v)
Origami ((1,2),(1,2))
```

"""
function origami(h::PermGroupElem, v::PermGroupElem)
  return Origami(GAP.Globals.Origami(GapObj(h), GapObj(v)), h, v)
end

function Base.show(io::IO, O::Origami)
  print(io, "Origami ($(O.h),$(O.v))")
end

GapObj(O::Origami) = O.o

function veech_group(O::Origami)
  GAP.Globals.VeechGroup(GapObj(O))
end