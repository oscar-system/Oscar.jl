global JTORIC_PATH = dirname(@__DIR__)

import GAP
import GAP: julia_to_gap, gap_to_julia

export GAP
export julia_to_gap, gap_to_julia

#using GAP

"""
    InstallDevelopmentVersionNConvex

Obtain a development version of NConvex.
"""
function InstallDevelopmentVersionNConvex()

    # find gap-location
    gap_location = GAP.GAPROOT
    
    # add pkg to this path
    p = splitpath( gap_location )
    p = replace.( p, "/"=>"")
    install_path = join( push!( p, "pkg" ), "/" )
    NConvex_path = join( push!( p, "NConvex" ), "/" )

    # inform what we are doing
    @info "Install development version of NConvex to \"" * NConvex_path * "\""
    
    # if this directory does exist, delete it
    if isdir( NConvex_path )

        # delete existing directory
        @info "DeleteExistingDirectory \"" * NConvex_path * "\""
        rm( NConvex_path, recursive=true )
        
    end
            
    # inform
    @info "Obtaining development version of NConvex into \"" * NConvex_path * "\""
        
    # prepare git operations
    res = GAP.Globals.LoadPackage(julia_to_gap("PackageManager"), false)
    @assert res
    git = julia_to_gap("git")
        
    # clone
    command = julia_to_gap( "clone https://github.com/HereAround/NConvex.git" )
    operation = GAP.Globals.PKGMAN_Exec(julia_to_gap( install_path ), git, command, julia_to_gap(""))
    if operation.code != 0
        @warn "Cloning of NConvex failed"
        return false
    end
    @info gap_to_julia(operation.output)
    
    # fetch development version
    command = julia_to_gap( "fetch origin martindevel:martindevel && git checkout martindevel" )
    operation = GAP.Globals.PKGMAN_Exec(julia_to_gap( NConvex_path ), git, command, julia_to_gap(""))
    if operation.code != 0
        @warn "Fetching of development version failed"
        return false
    end
    @info gap_to_julia(operation.output)
    
    # signal success
    return true
    
end
export InstallDevelopmentVersionNConvex
