###############################################################################
#
#  Temporary workarounds
#
###############################################################################
function symbols(Kt::Generic.RationalFunctionField)
    return Kt.S
end




###############################################################################
#
#  Includes
#
###############################################################################
include("semiring.jl")
include("semiring_map.jl")
include("matrix.jl")
include("poly.jl")
include("groebner_basis.jl")
include("initial.jl")
# include("groebner_polyhedron.jl")
# include("points.jl")
include("variety_supertype.jl")
include("hypersurface.jl")
include("curve.jl")
include("linear_space.jl")
include("variety.jl")
include("intersection.jl")
include("groebner_fan.jl")
