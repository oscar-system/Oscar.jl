module JToric


# use other Julia packages
using Pkg
using Oscar
using CapAndHomalg
using JToric

"""
Initializing function for 'JToric'.
"""
# initialization
function __init__()
    # define variable to call polymake
    JToric.GAP.Globals.POLYMAKE_JULIA_MODULE = Polymake
    
    # load necessary gap packages
    if ( ! GAP.Packages.load( "NConvex", "2021.04-24", install = false ) )
             warn("Could not load desired version of GAP package NConvex. JToric may not be fully functional.")
    end
    if ( ! GAP.Packages.load( "JConvex", "2021.08.12", install = false ) )
             warn("Could not load GAP package JConvex. JToric may not be fully functional.")
    end
    if ( ! GAP.Packages.load( "ToricVarieties", "2021.08.12", install = true ) )
             warn("Could not load nor install GAP package ToricVarieties. JToric may not be fully functional.")
    end
    
    # inform that JToric has been loaded
    print("Welcome to JToric ")
    printstyled("$version\n", color = :green)
    println("Enjoy this software to perform computations on toric geometry!")
    #println("Type: ?JToric for more information")
end


"""
    version

The version number of the loaded `JToric`.
Please mention this number in any bug report.
"""
global version = 0

if VERSION >= v"1.4"
    deps = Pkg.dependencies()
    if Base.haskey(deps, Base.UUID("9bfa4af0-0a41-40a3-86ae-7e8f6cc9e7ab"))
        ver = Pkg.dependencies()[Base.UUID("9bfa4af0-0a41-40a3-86ae-7e8f6cc9e7ab")]
        if occursin("/dev/", ver.source)
            version = VersionNumber("$(ver.version)-dev")
        else
            version = VersionNumber("$(ver.version)")
        end
    else
        version = "not installed"
    end
end

# include files
include("ToricVarieties.jl")
include("ToricDivisors.jl")
include("AttributesAndMethods.jl")

end
