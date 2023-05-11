export parametrization_plane_curve
export adjoint_ideal
export rational_point_conic
export parametrization_conic
export map_to_rational_normal_curve
export rat_normal_curve_anticanonical_map
export rat_normal_curve_It_Proj_Odd
export rat_normal_curve_It_Proj_Even
export invert_birational_map


################################################################################
function _tosingular(C::ProjPlaneCurve{QQFieldElem})
    F = C.eq
    T = parent(F)
    Tx = singular_poly_ring(T)
    return Tx(F)
end

function _fromsingular_ring(R::Singular.PolyRing)
    Kx = base_ring(R)
    if Kx isa Singular.N_AlgExtField
        FF, t = RationalFunctionField(QQ, "t")
        f = numerator(FF(Kx.minpoly))
        K, _ = number_field(f, "a")
    else
        K = QQ
    end
    newring, _ = polynomial_ring(K, symbols(R))
    return newring
end

function _tosingular_ideal(C::ProjCurve)
    I = C.I
    singular_assure(I)
    return I.gens.S
end

@doc raw"""
    parametrization_plane_curve(C::ProjPlaneCurve{QQFieldElem})

Return a rational parametrization of  `C`. 

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> C = ProjPlaneCurve(y^4-2*x^3*z+3*x^2*z^2-2*y^2*z^2)
Projective plane curve defined by -2*x^3*z + 3*x^2*z^2 + y^4 - 2*y^2*z^2

julia> parametrization_plane_curve(C)
3-element Vector{QQMPolyRingElem}:
 12*s^4 - 8*s^2*t^2 + t^4
 -12*s^3*t + 2*s*t^3
 8*s^4
```
"""
function parametrization_plane_curve(C::ProjPlaneCurve{QQFieldElem})
    s = "local"
    F = _tosingular(C)
    L = Singular.LibParaplanecurves.paraPlaneCurve(F, s)
    R = L[1]
    J = [L[2][i] for i in keys(L[2])][1]
    S = _fromsingular_ring(R)
    return gens(ideal(S, J))
end

@doc raw"""
    adjoint_ideal(C::ProjPlaneCurve{QQFieldElem})

Return the Gorenstein adjoint ideal of `C`. 

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> C = ProjPlaneCurve(y^4-2*x^3*z+3*x^2*z^2-2*y^2*z^2)
Projective plane curve defined by -2*x^3*z + 3*x^2*z^2 + y^4 - 2*y^2*z^2

julia> I = adjoint_ideal(C)
ideal(-x*z + y^2, x*y - y*z, x^2 - x*z)
```
"""
function adjoint_ideal(C::ProjPlaneCurve{QQFieldElem})
    n = 2
    F = _tosingular(C)
    R = parent(C.eq)
    I = Singular.LibParaplanecurves.adjointIdeal(F, n)
    return ideal(R, I)
end

@doc raw"""
    rational_point_conic(D::ProjPlaneCurve{QQFieldElem})

If the conic `D` contains a rational point, return the homogeneous coordinates of such a point.
If no such point exists, return a point on `D` defined over a quadratic field extension of $\mathbb Q$.
 
# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> D = ProjPlaneCurve(x^2 + 2*y^2 + 5*z^2 - 4*x*y + 3*x*z + 17*y*z);

julia> P = rational_point_conic(D)
3-element Vector{AbstractAlgebra.Generic.MPoly{nf_elem}}:
 -1//4*a
 -1//4*a + 1//4
 0

julia> S = parent(P[1])
Multivariate Polynomial Ring in x, y, z over Number field over Rational Field with defining polynomial t^2 - 2

julia> NF = base_ring(S)
Number field over Rational Field with defining polynomial t^2 - 2

julia> a = gen(NF)
a

julia> minpoly(a)
t^2 - 2
```
"""
function rational_point_conic(C::ProjPlaneCurve{QQFieldElem})
    F = _tosingular(C)
    L = Singular.LibParaplanecurves.rationalPointConic(F)
    R = L[1]
    P = [L[2][i] for i in keys(L[2])][1]
    S = _fromsingular_ring(R)
    return [S(P[1, i]) for i in 1:3]
end

@doc raw"""
    parametrization_conic(C::ProjPlaneCurve{QQFieldElem})

Given a conic `C`, return a vector `V` of polynomials in a new ring which should be
considered as the homogeneous coordinate ring of `PP^1`. The vector `V` defines a
rational parametrization `PP^1 --> C2 = {q=0}`.
"""
function parametrization_conic(C::ProjPlaneCurve{QQFieldElem})
    F = _tosingular(C)
    L = Singular.LibParaplanecurves.paraConic(F)
    R = L[1]
    J = [L[2][i] for i in keys(L[2])][1]
    S = _fromsingular_ring(R)
    return gens(ideal(S, J))
end

@doc raw"""
    map_to_rational_normal_curve(C::ProjPlaneCurve{QQFieldElem})

Return a rational normal curve of degree $\deg C-2$ which `C` is mapped.

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> C = ProjPlaneCurve(y^4-2*x^3*z+3*x^2*z^2-2*y^2*z^2);

julia> geometric_genus(C)
0

julia> map_to_rational_normal_curve(C)
Projective curve defined by the ideal(y(1)^2 + 2*y(1)*y(3) - 2*y(2)^2)
```
"""
function map_to_rational_normal_curve(C::ProjPlaneCurve{QQFieldElem})
    F = _tosingular(C)
    I = Singular.LibParaplanecurves.adjointIdeal(F)
    L = Singular.LibParaplanecurves.mapToRatNormCurve(F, I)
    S = _fromsingular_ring(L[1])
    J = [L[2][i] for i in keys(L[2])][1]
    IC = ideal(S, J)
    return ProjCurve(IC)
end


@doc raw"""
    rat_normal_curve_anticanonical_map(C::ProjCurve)

Return a vector `V` defining the anticanonical map `C --> PP^(n-2)`. Note that the
entries of `V` should be considered as representatives of elements in R/I,
where R is the basering.

# Examples
```jldoctest
julia> R, (v, w, x, y, z) = graded_polynomial_ring(QQ, ["v", "w", "x", "y", "z"])
(Multivariate Polynomial Ring in v, w, x, y, z over Rational Field graded by 
  v -> [1]
  w -> [1]
  x -> [1]
  y -> [1]
  z -> [1], MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[v, w, x, y, z])

julia> M = matrix(R, 2, 4, [v w x y; w x y z])
[v   w   x   y]
[w   x   y   z]

julia> V = minors(M, 2)
6-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 v*x - w^2
 v*y - w*x
 w*y - x^2
 v*z - w*y
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> RNC = ProjCurve(I)
Projective curve defined by the ideal(v*x - w^2, v*y - w*x, w*y - x^2, v*z - w*y, w*z - x*y, x*z - y^2)

julia> rat_normal_curve_anticanonical_map(RNC)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x
 -y
 z
```
"""
function rat_normal_curve_anticanonical_map(C::ProjCurve)
    R = base_ring(C.I)
    I = _tosingular_ideal(C)
    J = Singular.LibParaplanecurves.rncAntiCanonicalMap(I)
    return gens(ideal(R, J))
end

@doc raw"""
    rat_normal_curve_It_Proj_Odd(C::ProjCurve)

Return a vector `PHI` defining an isomorphic projection of `C` to `PP^1`.
Note that the entries of `PHI` should be considered as
representatives of elements in `R/I`, where `R` is the basering.

# Examples
```jldoctest
julia> R, (w, x, y, z) = graded_polynomial_ring(QQ, ["w", "x", "y", "z"]);

julia> M = matrix(R, 2, 3, [w x y; x y z])
[w   x   y]
[x   y   z]

julia> V = minors(M, 2)
3-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 w*y - x^2
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> TC = ProjCurve(I)
Projective curve defined by the ideal(w*y - x^2, w*z - x*y, x*z - y^2)

julia> rat_normal_curve_It_Proj_Odd(TC)
2-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 y
 -z
```
"""
function rat_normal_curve_It_Proj_Odd(C::ProjCurve)
    R = base_ring(C.I)
    I = _tosingular_ideal(C)
    J = Singular.LibParaplanecurves.rncItProjOdd(I)
    return gens(ideal(R, J))
end

# lookup an ideal with name s in the symbol table
# TODO move this to Singular.jl
function _lookup_ideal(R::Singular.PolyRingUnion, s::Symbol)
    for i in Singular.libSingular.get_ring_content(R.ptr)
        if i[2] == s
            @assert i[1] == Singular.mapping_types_reversed[:IDEAL_CMD]
            ptr = Singular.libSingular.IDEAL_CMD_CASTER(i[3])
            ptr = Singular.libSingular.id_Copy(ptr, R.ptr)
            return Singular.sideal{elem_type(R)}(R, ptr)
        end
    end
    error("could not find PHI")
end

@doc raw"""
    rat_normal_curve_It_Proj_Even(C::ProjCurve)

Return a vector `PHI` defining an isomorphic projection of `C` to `PP^1`.
Note that the entries of `PHI` should be considered as
representatives of elements in `R/I`, where `R` is the basering.

# Examples
```jldoctest
julia> R, (v, w, x, y, z) = graded_polynomial_ring(QQ, ["v", "w", "x", "y", "z"]);

julia> M = matrix(R, 2, 4, [v w x y; w x y z])
[v   w   x   y]
[w   x   y   z]

julia> V = minors(M, 2)
6-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 v*x - w^2
 v*y - w*x
 w*y - x^2
 v*z - w*y
 w*z - x*y
 x*z - y^2

julia> I = ideal(R, V);

julia> RNC = ProjCurve(I)
Projective curve defined by the ideal(v*x - w^2, v*y - w*x, w*y - x^2, v*z - w*y, w*z - x*y, x*z - y^2)

julia> rat_normal_curve_It_Proj_Even(RNC)
(MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, -y, z], Projective plane curve defined by -y(1)*y(3) + y(2)^2)
```
"""
function rat_normal_curve_It_Proj_Even(C::ProjCurve)
    R = base_ring(C.I)
    I = _tosingular_ideal(C)
    L = Singular.LibParaplanecurves.rncItProjEven(I)
    phi = _lookup_ideal(base_ring(I), :PHI)
    O = _fromsingular_ring(L[1]::Singular.PolyRing)
    conic = L[2][:CONIC]::Singular.spoly
    return gens(ideal(R, phi)), ProjPlaneCurve(O(conic))
end

@doc raw"""
    invert_birational_map(phi::Vector{T}, C::ProjPlaneCurve) where {T <: MPolyRingElem}

Return a dictionary where `image` represents the image of the birational map
given by `phi`, and `inverse` represents its inverse, where `phi` is a
birational map of the projective plane curve `C` to its image in the projective
space of dimension `size(phi) - 1`.
Note that the entries of `inverse` should be considered as
representatives of elements in `R/image`, where `R` is the basering.
"""
function invert_birational_map(phi::Vector{T}, C::ProjPlaneCurve) where {T <: MPolyRingElem}
    S = parent(phi[1])
    I = ideal(S, phi)
    singular_assure(I)
    L = Singular.LibParaplanecurves.invertBirMap(I.gens.S, _tosingular(C))
    R = _fromsingular_ring(L[1])
    J = L[2][:J]
    psi = L[2][:psi]
    return Dict([("image", gens(ideal(R, J))), ("inverse", gens(ideal(R, psi)))])
end
