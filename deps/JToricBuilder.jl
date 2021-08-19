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

    # add pkg to this path
    install_path = GAP.Packages.DEFAULT_PKGDIR
    NConvex_path = joinpath(install_path, "NConvex")

    # inform what we are doing
    @info "Install development version of NConvex to \"$(NConvex_path)\""

    # if this directory does exist, delete it
    if isdir( NConvex_path )
        # delete existing directory
        @info "DeleteExistingDirectory \"$(NConvex_path)\""
        rm( NConvex_path, recursive=true )
    end

    # inform
    @info "Obtaining development version of NConvex into \"$(NConvex_path)\""

    # clone
    run(`git clone -b martindevel https://github.com/HereAround/NConvex.git $NConvex_path`)

    # signal success
    return true
end
export InstallDevelopmentVersionNConvex
