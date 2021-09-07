export ProjCurve, defining_ideal

################################################################################
struct ProjCurve
    I::MPolyIdeal
    ambiant_dim::Int
    function ProjCurve(I::MPolyIdeal)
        dim(I) == 2 || error("wrong dimension for a projective curve")
        n = nvars(base_ring(I)) - 1
        new(I, n)
    end

    function Base.show(io::IO, C::ProjCurve)
        if !get(io, :compact, false)
            println(io, "Projective curve defined by the ", C.I)
        else
            println(io, C.I)
        end
    end
end

################################################################################

function Oscar.defining_ideal(C::ProjCurve)
    return C.I
end
