using Oscar
using Test

isdefined(Main, :test_save_load_roundtrip) ||
  include(joinpath(Oscar.oscardir, "test", "Serialization", "test_save_load_roundtrip.jl"))

include("session.jl")
include("basic_types.jl")
include("PolyhedralGeometry.jl")
include("containers.jl")
include("loading.jl")
include("ToricGeometry.jl")
include("Algebras.jl")
include("PolynomialsSeries.jl")
include("Matrices.jl")
include("Fields.jl")
include("TropicalGeometry.jl")
include("QuadForm.jl")
