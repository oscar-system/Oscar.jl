module JToric


# use other Julia packages
using Markdown
using Pkg
using Oscar

import Oscar:
    class_group,
    dim,
    iscomplete,
    isintegral,
    isnormal,
    isprincipal,
    issmooth,
    picard_group,
    projective_space


"""
Initializing function for 'JToric'.
"""
# initialization
function __init__()
    print("Welcome to JToric ")
    printstyled("$version\n", color = :green)
    println("Enjoy this software to perform computations on toric geometry!")
end

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])
const IS_DEV = isdir(@__DIR__, "..", ".git")

"""
    version

The version number of the loaded `JToric`.
Please mention this number in any bug report.
"""
const version = IS_DEV ? VersionNumber("$(VERSION_NUMBER)-dev") : VERSION_NUMBER

# include files
include("ToricVarieties.jl")
include("ToricDivisors.jl")
include("AttributesAndMethods.jl")

end
