@doc """
GroupLibraries2.jl provides access to several libraries of groups.
"""
module GroupLibraries2

# the necessary Julia packages
using Oscar
import Oscar.GAPGroup
import Oscar.IntegerUnion
import Oscar._permgroup_filter_attrs
import Oscar.translate_group_library_args
import Base.count, Base.get

# The following code can be loaded at compile time.
include("transitivegroups.jl")

export get, get_all, get_one, has, has_count, count, identification, has_identification
export TransitiveGroupsLibrary

end # module GroupLibraries2

using .GroupLibraries2

export get, get_all, get_one, has, has_count, count, identification, has_identification
export TransitiveGroupsLibrary
