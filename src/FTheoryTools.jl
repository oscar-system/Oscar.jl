module FTheoryTools

# use other Julia packages
using Markdown
using Pkg
using Oscar

"""
Initializing function for 'FTheoryTools'.
"""
# initialization
function __init__()
    print("Welcome to FTheoryTools ")
    printstyled("$version\n", color = :green)
    println("Enjoy this software to perform computations on resolutions of singular elliptic fibrations, with an eye towards F-theory.")
end

const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
const VERSION_NUMBER = VersionNumber(PROJECT_TOML["version"])
const IS_DEV = isdir(@__DIR__, "..", ".git")

"""
    version

The version number of the loaded `FTheoryTools`.
Please mention this number in any bug report.
"""
const version = IS_DEV ? VersionNumber("$(VERSION_NUMBER)-dev") : VERSION_NUMBER

# include files
include("WeierstrassModels/constructors.jl")
include("WeierstrassModels/attributes.jl")

include("TateModels/constructors.jl")
include("TateModels/attributes.jl")

include("TateModelsOverGeneralBaseSpace/constructors.jl")
include("TateModelsOverGeneralBaseSpace/attributes.jl")

end
