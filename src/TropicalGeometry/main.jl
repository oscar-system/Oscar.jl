_whyohwhy1 = AbstractAlgebra.Generic.MPolyRing{padic}
@attributes _whyohwhy1

_whyohwhy2 = FmpqMPolyRing
@attributes _whyohwhy2

_whyohwhy3 = FlintRationalField
@attributes _whyohwhy3


include("numbers.jl")
include("valuation.jl")
include("poly.jl")
include("initial.jl")


# Temporarily we will turn tropical polynomials into strings. This will be
# removed once Polymake.jl wraps its tropical polynomials and tropical numbers
#
# Warning: This function ignores all boundary cases!
function tropical_polynomial_to_polymake(f)
    convention = fun(base_ring(f))
    fstr = ""
    if convention == min
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


# Workaround for turning a PolyhedralFan of polymake into a proper PolyhedralComplex
function polyhedral_complex_workaround(pm::Polymake.BigObject)
    pc = pm
    typename = Polymake.type_name(pm)
    if typename[1:13] == "PolyhedralFan"
        pc = Polymake.fan.PolyhedralComplex(pm)
    end
    typename = Polymake.type_name(pc)
    if typename[1:17] != "PolyhedralComplex"
        error("Input object is not of type PolyhedralFan or PolyhedralComplex")
    end
    fv = Polymake.to_one_based_indexing(pc.FAR_VERTICES)
    mc = pc.MAXIMAL_POLYTOPES
    feasibles = [Polymake.to_zero_based_indexing(Polymake.row(mc, i)) for i in 1:Polymake.nrows(mc) if Polymake.incl(Polymake.row(mc, i), fv)>0]
    return Polymake.fan.PolyhedralComplex(POINTS=pc.VERTICES, INPUT_LINEALITY=pc.LINEALITY_SPACE, INPUT_POLYTOPES=feasibles)
end


include("variety_supertype.jl")
include("variety.jl")
include("hypersurface.jl")
include("curve.jl")
include("linear_space.jl")
