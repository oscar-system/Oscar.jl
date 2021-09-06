export parametrization_plane_curve, adjoint_ideal, rational_point_conic,
       parametrization_conic, map_to_rational_normal_curve,
       rat_normal_curve_anticanonical_map, rat_normal_curve_It_Proj_Odd,
       rat_normal_curve_It_Proj_Even


################################################################################
function _tosingular(C::ProjPlaneCurve{fmpq})
    F = C.eq
    T = parent(F)
    Tx = singular_ring(T)
    return Tx(F)
end

function _fromsingular_ring(R::Singular.PolyRing)
    Kx = base_ring(R)
    if typeof(Kx) == Singular.N_AlgExtField
        FF, t = RationalFunctionField(QQ, "t")
        f = numerator(FF(Kx.minpoly))
        K, _ = NumberField(f, "i")
    else
        K = QQ
    end
    newring, _ = PolynomialRing(K, [string(x) for x in gens(R)])
    return newring
end

function _tosingular_ideal(C::ProjCurve)
    I = C.I
    singular_assure(I)
    return I.gens.S
end

@doc Markdown.doc"""
    parametrization_plane_curve(C::ProjPlaneCurve{fmpq}, s::String = "local")

Return a vector `V` of polynomials in a new ring `R` consisting of a rational parametrization
of the rational plane curve `C`. The entries of `V` parametrize the three
coordinates of the rational curve. The ground field of `R` is either QQ or some
algebraic extension QQ(a). The optional string can be "normal", to compute the
integral basis via normalization, or local (default), to make local analysis of
singularities first and apply normalization separately.
"""
function parametrization_plane_curve(C::ProjPlaneCurve{fmpq}, s::String = "local")
    F = _tosingular(C)
    L = Singular.LibParaplanecurves.paraPlaneCurve(F, s)
    R = L[1]
    J = [L[2][i] for i in keys(L[2])][1]
    S = _fromsingular_ring(R)
    return gens(ideal(S, J))
end

@doc Markdown.doc"""
    adjoint_ideal(C::ProjPlaneCurve{fmpq}, n::Int = 2)

Return the adjoint ideal of the curve `C`. If n = 1, the computation is done via
normalization, if n = 2 (default value), the local analysis of singularities is
made first and then the normalization is applied separately, if n = 3 (resp
n = 4), it uses normalization via (resp local) ideal quotient.
"""
function adjoint_ideal(C::ProjPlaneCurve{fmpq}, n::Int = 2)
    F = _tosingular(C)
    R = parent(C.eq)
    I = Singular.LibParaplanecurves.adjointIdeal(F, n)
    return ideal(R, I)
end

@doc Markdown.doc"""
    rational_point_conic(C::ProjPlaneCurve{fmpq})

Given a conic `C`, return the coordinates of a point on it. If there is no
rational point on the curve `C`, it creates a field extension of `QQ`.
"""
function rational_point_conic(C::ProjPlaneCurve{fmpq})
    F = _tosingular(C)
    L = Singular.LibParaplanecurves.rationalPointConic(F)
    R = L[1]
    P = [L[2][i] for i in keys(L[2])][1]
    S = _fromsingular_ring(R)
    return [S(P[1, i]) for i in 1:3]
end

@doc Markdown.doc"""
    parametrization_conic(C::ProjPlaneCurve{fmpq})

Given a conic `C`, return a vector `V` of polynomials in a new ring which should be
considered as the homogeneous coordinate ring of `PP^1`. The vector `V` defines a
rational parametrization `PP^1 --> C2 = {q=0}`.
"""
function parametrization_conic(C::ProjPlaneCurve{fmpq})
    F = _tosingular(C)
    L = Singular.LibParaplanecurves.paraConic(F)
    R = L[1]
    J = [L[2][i] for i in keys(L[2])][1]
    S = _fromsingular_ring(R)
    return gens(ideal(S, J))
end

@doc Markdown.doc"""
    map_to_rational_normal_curve(C::ProjPlaneCurve)

Return a rational normal curve to which the curve `C` is mapped.
"""
function map_to_rational_normal_curve(C::ProjPlaneCurve{fmpq})
    F = _tosingular(C)
    I = Singular.LibParaplanecurves.adjointIdeal(F)
    L = Singular.LibParaplanecurves.mapToRatNormCurve(F, I)
    S = _fromsingular_ring(L[1])
    J = [L[2][i] for i in keys(L[2])][1]
    IC = ideal(S, J)
    return ProjCurve(IC)
end


@doc Markdown.doc"""
    rat_normal_curve_anticanonical_map(C::ProjCurve)

Return a vector `V` defining the anticanonical map `C --> PP^(n-2)`. Note that the
entries of `V` should be considered as representatives of elements in R/I,
where R is the basering.
"""
function rat_normal_curve_anticanonical_map(C::ProjCurve)
    R = base_ring(C.I)
    I = _tosingular_ideal(C)
    J = Singular.LibParaplanecurves.rncAntiCanonicalMap(I)
    return gens(ideal(R, J))
end

@doc Markdown.doc"""
    rat_normal_curve_It_Proj_Odd(C::ProjCurve)

Return a vector `PHI` defining an isomorphic projection of `C` to `PP^1`.
Note that the entries of `PHI` should be considered as
representatives of elements in `R/I`, where `R` is the basering.
"""
function rat_normal_curve_It_Proj_Odd(C::ProjCurve)
    R = base_ring(C.I)
    I = _tosingular_ideal(C)
    J = Singular.LibParaplanecurves.rncItProjOdd(I)
    return gens(ideal(R, J))
end

function rat_normal_curve_It_Proj_Even(C::ProjCurve)
    I = _tosingular_ideal(C)
    L = Singular.LibParaplanecurves.rncItProjEven(I)
    S = _fromsingular_ring(L[1])
    f = [L[2][i] for i in keys(L[2])][1]
    return ProjPlaneCurve(S(f))
end
