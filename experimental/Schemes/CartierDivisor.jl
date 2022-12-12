export CartierDivisor

@attributes mutable struct CartierDivisor{
                                          CoveredSchemeType<:AbsCoveredScheme
                                         }
    X::CoveredSchemeType
    I::IdealSheaf

    function CartierDivisor(
            X::AbsCoveredScheme, 
            D::IdDict{<:AbsSpec, <:RingElem};
            check::Bool=true
        )
        all(x->(x in keys(D)), affine_charts(X)) || error("equations must be given on all affine charts of the variety")
        ID = IdDict{AbsSpec, Ideal}()
        for U in keys(D)
            ID[U] = ideal(OO(U), D[U])
        end
        I = IdealSheaf(X, ID)
        if check
            # TODO: Check the cocycle condition
        end
        return new{typeof(X)}(X, I)
    end
end

function (C::CartierDivisor)(U::AbsSpec)
    return gens(C.I(U))[1]
end

variety(C::CartierDivisor) = C.X
scheme(C::CartierDivisor) = C.X
