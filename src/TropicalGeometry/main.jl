
include("semiring.jl")
include("valuation.jl")
include("poly.jl")
include("initial.jl")
include("groebner_basis.jl")
include("groebner_polyhedron.jl")
include("points.jl")

# Temporarily we will turn tropical polynomials into strings. This will be
# removed once Polymake.jl wraps its tropical polynomials and tropical numbers
#
# Warning: This function ignores all boundary cases!
function tropical_polynomial_to_polymake(f)
    fstr = ""
    if convention(base_ring(f)) == min
        fstr *= "min("
    else
        fstr *= "max("
    end
    td = total_degree(f)
    for i in 1:length(f)
        fstr *= repr(coeff(f,i).data)
        e = exponent_vector(f,i)
        if td - sum(e) != 0
            fstr *= "+"
            fstr *= repr(td-sum(e))
            fstr *= "x_0"
        end
        if !iszero(e)
            for j in 1:length(e)
                if !iszero(e[j])
                    fstr *= "+"
                    fstr *= repr(e[j])
                    fstr *= "*x_"
                    fstr *= repr(j)
                end
            end
        end
        if i != length(f)
            fstr *= ","
        end
    end
    fstr *= ")"
    result = ["x_"*repr(i) for i in 0:nvars(parent(f))]
    prepend!(result, [fstr])
    return result
end


###
# Allow gcd of vectors of univariate rational polynomials
# to make their handling similar to that of integers
###
gcd(F::Vector{QQPolyRingElem}) = reduce(gcd, F)
gcd(F::QQPolyRingElem...) = reduce(gcd, F)



include("variety_supertype.jl")
include("variety.jl")
include("hypersurface.jl")
include("curve.jl")
include("linear_space.jl")
