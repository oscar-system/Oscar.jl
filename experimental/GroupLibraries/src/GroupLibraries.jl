@doc """
GroupLibraries.jl provides access to several libraries of groups.
"""
module GroupLibraries

# the necessary Julia packages
using Oscar

# The following code can be loaded at compile time.
include("transitivegroups.jl")

export TransitiveGroups

end # module GroupLibraries

using .GroupLibraries

export TransitiveGroups
