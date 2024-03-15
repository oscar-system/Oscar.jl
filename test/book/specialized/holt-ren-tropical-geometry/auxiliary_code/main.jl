using Pkg
Pkg.add("MixedSubdivisions"; io=devnull)

import LibGit2
LibGit2.clone("https://github.com/isaacholt100/generic_root_count", "generic_root_count")
include(joinpath(@__DIR__,"generic_root_count","src","main.jl"));
