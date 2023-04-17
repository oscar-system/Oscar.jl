using Oscar
using Test

include("gap_to_oscar.jl")
include("iso_gap_oscar.jl")
include("iso_oscar_gap.jl")
include("oscar_to_gap.jl")

# Test the OscarInterface package
@test GAP.Globals.TestDirectory(GAP.Globals.DirectoriesPackageLibrary(GAP.Obj("OscarInterface"), GAP.Obj("tst")))
