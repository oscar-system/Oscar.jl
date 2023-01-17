export CartierDivisor
export trivializing_covering

@attributes mutable struct CartierDivisor{
                                          CoveredSchemeType<:AbsCoveredScheme
                                         }
    X::CoveredSchemeType
    I::IdealSheaf
    C::Covering

    function CartierDivisor(
            X::AbsCoveredScheme, 
            D::IdDict{<:AbsSpec, <:RingElem};
            check::Bool=true
        )
        ID = IdDict{AbsSpec, Ideal}()
        for U in keys(D)
            ID[U] = ideal(OO(U), D[U])
        end
        I = IdealSheaf(X, ID, check=check)
        if check
          for U in keys(D)
            is_zero_divisor(OO(U)(D[U])) && error("local elements must not be zero divisors")
          end
          # TODO: 
          # - Check that every affine chart is covered
        end
        C = Covering(collect(keys(D)), IdDict{Tuple{AbsSpec, AbsSpec}, AbsGlueing}())
        return new{typeof(X)}(X, I, C)
    end
end

function (C::CartierDivisor)(U::AbsSpec)
    return gens(C.I(U))
end

variety(C::CartierDivisor) = C.X
scheme(C::CartierDivisor) = C.X
trivializing_covering(C::CartierDivisor) = C.C

