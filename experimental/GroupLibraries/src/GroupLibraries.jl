@doc """
GroupLibraries.jl provides access to several libraries of groups.
"""
module GroupLibraries

# the necessary Julia packages
using Oscar

# The following code can be loaded at compile time.
include("transitivegroups.jl")
include("perfectgroups.jl")

export TransitiveGroups
export PerfectGroups

end # module GroupLibraries

using .GroupLibraries

export TransitiveGroups
export PerfectGroups
