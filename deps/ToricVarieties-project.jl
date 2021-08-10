global JTORIC_PATH = dirname(@__DIR__)

import CapAndHomalg
import GAP
import GAP: julia_to_gap, gap_to_julia, @g_str, @gap, GapObj

export CapAndHomalg
export GAP
export julia_to_gap, gap_to_julia, @g_str, @gap, GapObj

"""
    JToric.CAPANDHOMALG_PKG_DIR

The directory in which CapAndHomalg install [GAP](https://www.gap-system.org)
packages from the GitHub organization
[https://github.com/homalg-project/](https://github.com/homalg-project/).
"""
function CapAndHomalgPath()

    p = splitpath( pathof(CapAndHomalg) )
    p = replace.(p, "/"=>"")
    pop!( p )
    pop!( p )
    push!( p , "pkg" )    
    return join( p, "/" )

end
global CAPANDHOMALG_PKG_DIR = CapAndHomalgPath()


"""
    Build_ToricVarieties_project

Builds the software of the ToricVarieties_project.
"""
function BuildToricVarieties_project()

        # prepare
        script = julia_to_gap("./build.sh")
        dir = joinpath(CAPANDHOMALG_PKG_DIR, "ToricVarieties_project")
        if !isdir(dir)
            return false
        end

        # build ToricVarieties_project
        @info "Building ToricVarieties"
        operation = GAP.Globals.PKGMAN_Exec(julia_to_gap(dir), script, julia_to_gap(""), julia_to_gap(""))
        if operation.code != 0
            @warn "Building failed:\n" * gap_to_julia(operation.output)
            return false
        end
        @info gap_to_julia(operation.output)
                
        # signal success
        print( "\n\n" )
        println("Done building ToricVarieties_project.")
        return true;

end
export BuildToricVarieties_project


"""
    UpdateNConvex

Updates NConvex to a version compatible with JConvex and JToric.
"""
function UpdateNConvex()

        # set up git and the pull command
        res = GAP.Globals.LoadPackage(julia_to_gap("PackageManager"), false)
        @assert res
        git = julia_to_gap("git")
        command = julia_to_gap( "fetch git://github.com/HereAround/NConvex.git martindevel && git merge FETCH_HEAD" )
        
        # save directory of NConvex
        dir = joinpath(CAPANDHOMALG_PKG_DIR, "NConvex")

        # check if this directory does exist
        if !isdir(dir)
            return false
        end

        # inform that we update NConvex now
        @info "Updating \"" * dir * "\""
        
        # add remote
        operation = GAP.Globals.PKGMAN_Exec(julia_to_gap(dir), git, command, julia_to_gap(""))
        if operation.code != 0
            @warn "Fetching development version failed:\n" * gap_to_julia(operation.output)
            return false
        end
        @info gap_to_julia(operation.output)
        
        # signal success
        return true

end
export UpdateNConvex
