module JToric


# use other Julia packages
using Markdown
using Pkg
using Oscar
using CapAndHomalg

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
    # expose the Polymake module to GAP for use in JConvex
    GAP.Globals._Polymake_jl = Oscar.Polymake

    # ensure "our" JConvex is loaded
    GAP.Globals.SetPackagePath(GapObj("JConvex"), GapObj(abspath(@__DIR__, "..", "pkg", "JConvex")))

    # load necessary gap packages
    if ( ! GAP.Packages.load( "NConvex", "2021.04-24", install = false ) )
             @warn("Could not load desired version of GAP package NConvex. JToric may not be fully functional.")
    end
    if ( ! GAP.Packages.load( "JConvex", "2021.09.21", install = false ) )
             @warn("Could not load GAP package JConvex. JToric may not be fully functional.")
    end
    if ( ! GAP.Packages.load( "ToricVarieties", "2021.08.12", install = true ) )
             @warn("Could not load nor install GAP package ToricVarieties. JToric may not be fully functional.")
    end
    
    # inform that JToric has been loaded
    print("Welcome to JToric ")
    printstyled("$version\n", color = :green)
    println("Enjoy this software to perform computations on toric geometry!")
    #println("Type: ?JToric for more information")
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
include("conversion.jl")
include("ToricVarieties.jl")
include("ToricDivisors.jl")
include("AttributesAndMethods.jl")

end
