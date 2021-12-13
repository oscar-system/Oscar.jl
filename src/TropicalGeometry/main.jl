include("Tropical.jl")

# Temporarily we will turn tropical polynomials into strings. This will be
# removed once Polymake.jl wraps its tropical polynomials and tropical numbers
#
# Warning: This function ignores all boundary cases!
function tropical_polynomial_to_polymake(f)
    convention = fun(base_ring(f))
    result = ""
    if convention == min
        result *= "min("
    else
        result *= "max("
    end
    for i in 1:length(f)
        result *= repr(coeff(f,i).data) 
        e = exponent_vector(f,i)
        if !iszero(e)
            for j in 1:length(e)
                if !iszero(e[j])
                    result *= "+"
                    result *= repr(e[j])
                    result *= "*x("
                    result *= repr(j)
                    result *= ")"
                end
            end
        end
        if i != length(f)
            result *= ","
        end
    end
    result *= ")"
end

include("TropicalVarieties.jl")
