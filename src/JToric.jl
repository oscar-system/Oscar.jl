module JToric


# load GAP and CapAndHomalg
import CapAndHomalg
import GAP

# use packages
using Pkg
using JToric

"""
Initializing function for 'JToric'.
"""
# initialization
function __init__()
    # load necessary gap packages
    CapAndHomalg.LoadPackage( "JuliaInterface" )
    CapAndHomalg.LoadPackage( "JConvex" )
    CapAndHomalg.LoadPackage( "ToricV" )
    
    # inform that JToric has been loaded
    print("Welcome to JToric ")
    printstyled("$version\n", color = :green)
    println("Enjoy this software to perform computations on toric geometry!")
    #println("Type: ?JToric for more information")
end


"""
    JToric.version

The version number of the loaded `CapAndHomalg`.
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
include("AttributesAndMethods.jl")

end
