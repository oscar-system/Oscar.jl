###############################################################################
#
#  Temporary workarounds
#
###############################################################################
function symbols(Kt::Generic.RationalFunctionField)
    return Kt.S
end

export Min, Max

"""
    Min

Indicating that we are in the tropical semiring ``(ℝ ∪ {∞}, ⊕, ⊙)`` with the
min convention.
"""
struct Min end

"""
    Max

Indicating that we are in the tropical semiring ``(ℝ ∪ {-∞}, ⊕, ⊙)`` with the
max convention.
"""
struct Max end

const MinMax = Union{Min, Max}

# to make Max(x, y) and Min(x, y) work
(::Type{Min})(x, y) = min(x, y)

(::Type{Max})(x, y) = max(x, y)

convention_function(::Type{Min}) = min

convention_function(::Type{Max}) = max



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
