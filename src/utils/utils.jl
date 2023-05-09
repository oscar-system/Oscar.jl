const cornerstones = String["AbstractAlgebra", "GAP", "Hecke", "Nemo", "Polymake", "Singular"];
const jll_deps = String["Antic_jll", "Arb_jll", "Calcium_jll", "FLINT_jll", "GAP_jll",
                        "libpolymake_julia_jll", "libsingular_julia_jll",
                        "polymake_jll", "Singular_jll"];

include("versioninfo.jl")
include("docs.jl")
include("tests.jl")
include("rand.jl")
