export sum_Point_EllCurveZnZ,
    ECM,
    rand_pair_EllCurve_Point,
    IntMult_Point_EllCurveZnZ,
    cornacchia_algorithm,
    Miller_Rabin_test,
    Pollard_rho,
    Pollard_p_1,
    ECPP

################################################################################
# Elliptic curves over a ring Z/nZ
################################################################################
################################################################################
# Helping functions
################################################################################
# addition on two projective points on an elliptic curve over Z/nZ, given in
# short Weierstrass form, when it is possible. It gives as an output a pair
# composed of the coordinates of the sum and ZZ(1) if the sum exists, and
# (0, 0, 0) and the gcd otherwise.
# Adapted from the corresponding code in Hecke.

function _add(
    P::Vector{Nemo.fmpz_mod},
    Q::Vector{Nemo.fmpz_mod},
    E::Vector{Nemo.fmpz_mod},
)
    length(P) == 3 && length(Q) == 3 || error("arrays of size 3 required")
    length(E) == 2 || error("array of size 2 required")
    A = parent(P[1])
    all(x -> parent(x) == A, [P..., Q..., E...]) || error("Not the same parent")
    P[3] != A(1) &&
        P != [A(0), A(1), A(0)] &&
        error("require infinity point or last coordinate 1")
    Q[3] != A(1) &&
        Q != [A(0), A(1), A(0)] &&
        error("require infinity point or last coordinate 1")
    n = modulus(A)
    if P == [A(0), A(1), A(0)]
        return [Q, ZZ(1)]
    elseif Q == [A(0), A(1), A(0)]
        return [P, ZZ(1)]
    end
    if P[1] != Q[1]
        d = Q[1] - P[1]
        g = gcd(n, Hecke.data(d))
        if g == ZZ(1)
            m = divexact(Q[2] - P[2], d)
            x = m^2 - P[1] - Q[1]
            y = m * (P[1] - x) - P[2]
            z = A(1)
        else
            x = A(0)
            y = A(0)
            z = A(0)
        end
    elseif P[2] != Q[2]
        return [[A(0), A(1), A(0)], ZZ(1)]
    elseif P[2] != 0
        d = 2 * P[2]
        g = gcd(n, Hecke.data(d))
        if g == ZZ(1)
            m = divexact(3 * (P[1])^2 + E[1], 2 * P[2])
            x = m^2 - 2 * P[1]
            y = m * (P[1] - x) - P[2]
            z = A(1)
        else
            x = A(0)
            y = A(0)
            z = A(0)
        end
    else
        return [[A(0), A(1), A(0)], ZZ(1)]
    end
    return [[x, y, z], g]
end

################################################################################
# Creates a pair composed of the coordinates of a point and the coefficients of
# an elliptic curve going through that point.

function _rand_point_curve(A::Nemo.FmpzModRing)
    n = modulus(A)
    a = A(rand(ZZ(1):n))
    P = [A(rand(ZZ(1):n)), A(rand(ZZ(1):n)), A(ZZ(1))]
    b = A(P[2]^2 - P[1]^3 - a * P[1])
    return [P, [a, b]]
end

################################################################################
# Creates a list of pairs composed of a the coefficitents of an elliptic curve
# and the coordinates of a point on it.

function _rand_list(A::Nemo.FmpzModRing, N::Int)
    E = Vector{Nemo.fmpz_mod}[]
    P = Vector{Nemo.fmpz_mod}[]
    for i = 1:N
        L = _rand_point_curve(A)
        E = push!(E, L[2])
        P = push!(P, L[1])
    end
    return [E, P]
end

################################################################################
# Return the coordinates of m*P and ZZ(1) when the computation is possible, and
# returns (0, 0, 0) and the gcd otherwise.

function _scalar_mult(P::Vector{Nemo.fmpz_mod}, E::Vector{Nemo.fmpz_mod}, m::fmpz)
    length(P) == 3 || error("arrays of size 3 required")
    length(E) == 2 || error("array of size 2 required")
    A = parent(P[1])
    res = [[A(0), A(1), A(0)], ZZ(1)]
    if m < 0
        error("Positive integer expected")
    elseif m == 0
        return res
    elseif m == 1
        return [P, ZZ(1)]
    else
        b = m
        Q = [P, ZZ(1)]
        while !iszero(b)
            if b % 2 == 1
                res = _add(res[1], Q[1], E)
                res[2] != ZZ(1) && return res
            end
            Q = _add(Q[1], Q[1], E)
            Q[2] != ZZ(1) && return Q
            b = div(b, 2)
        end
        return res
    end
end

################################################################################
# Return the elliptic plane curve corresponding to the array.

function _toProjEllipticCurve(R::MPolyRing{S}, E::Vector{S}) where {S<:Nemo.fmpz_mod}
    length(E) == 2 || error("array of size 2 required")
    x = gen(R, 1)
    y = gen(R, 2)
    z = gen(R, 3)
    return ProjEllipticCurve(R(y^2 * z - x^3 - E[1] * x * z^2 - E[2] * z^3))
end

################################################################################
# Functions
################################################################################
@doc Markdown.doc"""
    ECM(n::fmpz; nbcurve::Int = 25000, multfact::fmpz = factorial(ZZ(10^4)))

Return a factor of `n`, obtained with the Elliptic Curve Method.
"""
function ECM(n::fmpz; nbcurve::Int = 25000, multfact::fmpz = factorial(ZZ(10^4)))
    A = ResidueRing(ZZ, n)
    N = nbcurve
    L = _rand_list(A, N)

    for i = 1:length(L)
        Q = _scalar_mult(L[2][i], L[1][i], multfact)
        Q[2] > 1 && return Q[2]
    end
    return ZZ(1)
end

################################################################################
# the sum of two points might not be a point, this is not a group operation.

@doc Markdown.doc"""
    sum_Point_EllCurveZnZ(P::Point_EllCurve{S}, Q::Point_EllCurve{S}) where S <: Nemo.fmpz_mod

Return, if possible, the sum of the points `P` and `Q`, and an error otherwise.
"""
function sum_Point_EllCurveZnZ(
    P::Point_EllCurve{S},
    Q::Point_EllCurve{S},
) where {S<:Nemo.fmpz_mod}
    A = parent(P.Pt[1])
    E = P.C
    E.Hecke_ec.short || error("requires short Weierstrass form")
    Q = _add(P.Pt.v, Q.Pt.v, E.Hecke_ec.coeff)
    PP = Oscar.Geometry.parent(P.Pt)
    if Q[2] == A(1)
        return Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, Q[1]))
    else
        error(Q[2], " is not invertible in the base ring, cannot perform the sum")
    end
end

################################################################################
@doc Markdown.doc"""
    IntMult_Point_EllCurveZnZ(m::fmpz, P::Point_EllCurve{S}) where S <: Nemo.fmpz_mod

Return, if possible, the point `mP`, and an error otherwise.
"""
function IntMult_Point_EllCurveZnZ(m::fmpz, P::Point_EllCurve{S}) where {S<:Nemo.fmpz_mod}
    E = P.C
    E.Hecke_ec.short || error("requires short Weierstrass form")
    PP = Oscar.Geometry.parent(P.Pt)
    Q = _scalar_mult(P.Pt.v, E.Hecke_ec.coeff, m)
    if Q[2] == ZZ(1)
        return Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, Q[1]))
    else
        error(Q[2], " is not invertible in the base ring, cannot perform the sum")
    end
end

################################################################################
@doc Markdown.doc"""
    rand_pair_EllCurve_Point(R::Oscar.MPolyRing_dec{S}, PP::Oscar.Geometry.ProjSpc{S}) where S <: Nemo.fmpz_mod

Return a pair composed of an elliptic plane curve `E` with equation in `R`,
and a point `P` on `E`.
"""
function rand_pair_EllCurve_Point(
    R::Oscar.MPolyRing_dec{S},
    PP::Oscar.Geometry.ProjSpc{S},
) where {S<:Nemo.fmpz_mod}
    A = base_ring(R)
    n = modulus(A)
    i = 0
    L = _rand_point_curve(A)
    while i < 100 && !is_unit(4 * L[2][1]^3 + 27 * L[2][2]^2)
        i = i + 1
        L = _rand_list(A, 1)
    end
    if is_unit(4 * L[2][1]^3 + 27 * L[2][2]^2)
        E = _toProjEllipticCurve(R, L[2])
        return [E, Point_EllCurve(E, Oscar.Geometry.ProjSpcElem(PP, L[1]))]
    else
        error("Did not manage do find an elliptic curve with invertible discriminant")
    end
end

################################################################################
# ECPP
################################################################################
# Cornacchia Algorithm
################################################################################

function evenodd(n::fmpz, fac::fmpz = ZZ(2))
    res = (ZZ(1), n)
    k = fac
    g = gcd(k, n)
    while g != ZZ(1)
        if g < k
            return (res[1] * g, div(res[2], g))
        end
        res = (res[1] * k, div(res[2], k))
        k = k^2
        g = gcd(k, res[2])
    end
    return res
end

################################################################################
@doc Markdown.doc"""
    cornacchia_algorithm(d::fmpz, m::fmpz)

Return `true` and a solution of `x^2 + d*y^2 = m` if it exists, and false and
`(0, 0)` otherwise.
"""
function cornacchia_algorithm(d::fmpz, m::fmpz)
    (even, odd) = evenodd(m)
    D = abs(d)
    R = ResidueRing(ZZ, odd)
    S = ResidueRing(ZZ, min(ZZ(8), even))
    X = issquare_with_square_root(-R(D))
    if !X[1]
        return (false, (ZZ(0), ZZ(0)))
    end
    s1 = Hecke.data(X[2])
    if !isone(even)
        !isone(S(D)) && (false, (ZZ(0), ZZ(0)))
        s2 = sqrtmod(-D, even)
        s = abs(crt(-s1, odd, -s2 + ZZ(2)^(Int(log2(even) - 1)), even, true))
    else
        s = s1
    end

    if 2 * s > m
        B = m - s
    else
        B = s
    end

    A = m
    while B^2 >= m
        A, B = B, A % B
    end
    t = m - B^2
    L = divrem(t, D)

    if !iszero(L[2]) || !is_square(L[1])
        return (false, (ZZ(0), ZZ(0)))
    end

    return (true, [B, isqrt(L[1])])
end

################################################################################
###############Hilbert Class Polynomial#########################################
################################################################################

# Compute reduced forms of imaginary quadratic number fields
# based on Algo 5.3.5 in the book by Cohen
function class_number(D::fmpz)
    D < ZZ(0) || error("Input needs to be negative")

    h = 0
    b = abs(D % 2)
    B = root(div(-D, 3), 2)

    while b <= B
        q = div(b^2 - D, 4)
        a = ZZ(max(b, 1))
        while a^2 <= q
            if q % a == 0 && gcd(a, b, div(q, a)) == 1
                if a == b || a^2 == q || b == 0
                    h = h + 1
                else
                    h = h + 2
                end
            end
            a = a + 1
        end
        b = b + 2
    end
    return ZZ(h)
end

#################################################################################

function funddiscriminant(n::Int)
    D = ZZ(-3)
    L = Vector{fmpz}()
    while abs(D) <= n
        if mod(D, 4) == 1 && is_squarefree(D)
            h = class_number(D)
            push!(L, D)
        end
        if mod(D, 4) == 0 &&
           is_squarefree(div(D, 4)) &&
           (mod(div(D, 4), 4) == 2 || mod(div(D, 4), 4) == 3)
            h = class_number(D)
            push!(L, D)
        end
        D = D - 1
    end

    return L
end

################################################################################
# Miller Rabin Test
################################################################################
function _Miller_Rabin_witness(N::fmpz, a::fmpz, even::fmpz, odd::fmpz)
    R = ResidueRing(ZZ, N)
    x = ZZ(R(a)^odd)
    if x == ZZ(1) || x == N - ZZ(1)
        return false
    else
        for i = 1:nbits(even)
            x = ZZ(R(x)^2)
            if x == N - ZZ(1)
                return false
            end
        end
    end
    return true
end

################################################################################

@doc Markdown.doc"""
    Miller_Rabin_test(N::fmpz, k::Int64 = 20)

Given an odd number `N`, return `false` if the number is composite, and `true` if it
is probably prime.
"""
function Miller_Rabin_test(N::fmpz, k::Int64 = 20)
    !iszero(N % 2) || error("an odd number is expected")
    (even, odd) = evenodd(N - ZZ(1))
    for i = 1:k
        a = rand(ZZ(2):min(N, ZZ(10^6)))
        if _Miller_Rabin_witness(N, a, even, odd)
            return false
        end
    end
    return true
end

################################################################################
# Pollard's methods
################################################################################
@doc Markdown.doc"""
    Pollard_rho(N::fmpz, bound::Int = 50000)

The algorithm computes a factor of `N` using the Pollard rho algorithm
and returns it.
"""
function Pollard_rho(N::fmpz, bound::Int = 50000)
    R = ResidueRing(ZZ, N)
    x = rand(ZZ(2):N-1)
    y = rand(ZZ(2):N-1)
    d = ZZ(1)
    i = 1
    while d == 1 && i <= bound
        x = ZZ(R(x)^2 + R(1))
        y = R(y)^2 + R(1)
        y = y + R(1)
        d = gcd(x - ZZ(y), N)
        i = i + 1
    end
    return d
end

################################################################################

@doc Markdown.doc"""
    Pollard_p_1(N::fmpz, B::fmpz = ZZ(10)^5)

The algorithm computes a factor of `N` and returns it.
"""
function Pollard_p_1(N::fmpz, B::fmpz = ZZ(10)^5)
    p = [i for i in PrimesSet(ZZ(1), B)]
    x = rand(ZZ(2):N-1)
    y = x
    c = ZZ(0)
    i = 0
    j = i
    k = length(p)

    while true
        i = i + 1
        if i > k
            g = gcd(x - ZZ(1), N)
            if g == ZZ(1)
                return g
            else
                i = j
                x = y
                break
            end
        else
            q = p[i]
            qq = q
            l = div(B, q)
            while qq <= ZZ(l)
                qq = q * qq
            end
            x = x^qq % N
            c = c + ZZ(1)
            if c >= ZZ(20)
                g = gcd(x - ZZ(1), N)
                if g == ZZ(1)
                    c = ZZ(0)
                    j = i
                    y = x
                else
                    i = j
                    x = y
                    break
                end
            end
        end
    end

    while i <= k

        i = i + 1
        q = p[i]
        qq = q
        t = 0
        while true

            x = x^q % N
            g = gcd(x - ZZ(1), N)
            !isone(g) && return g
            qq = q * qq
            if qq > B
                break
            end
        end
    end
    return gcd(x - ZZ(1), N)
end

################################################################################
############ ECPP helper functions #############################################
################################################################################

function trial_division(N::fmpz)
    D = Hecke.factor_trial_range(N)[1]
    V = [p^D[p] for p in keys(D)]
    t = prod(v for v in V)
    m = div(N, t)
    m != 1 && return (t, m)
    return (ZZ(1), maximum(keys(D)))
end

################################################################################

function find_q(N::fmpz, D::fmpz)
    (bool, (x, y)) = cornacchia_algorithm(D, N)
    if !bool
        return (false, (ZZ(1), ZZ(N)))
    end

    if D == ZZ(-3)
        M = [
            N + 1 + 2 * x,
            N + 1 - 2 * x,
            N + 1 + x + 3 * y,
            N + 1 - x - 3 * y,
            N + 1 - x + 3 * y,
            N + 1 + x - 3 * y,
        ]
    elseif D == ZZ(-4)
        M = [N + 1 + 2 * x, N + 1 - 2 * x, N + 1 + 4 * y, N + 1 - 4 * y]
    else
        M = [N + 1 + 2 * x, N + 1 - 2 * x]
    end
    for m in M
        # split even and odd part
        (_, q) = evenodd(m)

        #trial divisions
        (_, q) = trial_division(q)

        #Factoring via ECM
        t = ZZ(2)
        i = 0

        while t != ZZ(1) && !Miller_Rabin_test(q) && i <= 10
            i = i + 1
            t = ZZ(1)
            t = lcm(t, ECM(q))
            if t == ZZ(1)
                t = ZZ(2)
            else
                (_, q) = evenodd(q, t)
            end
        end
        sqrt(BigInt(q)) > sqrt(sqrt(BigInt(N))) + 1 && return (true, (m, q))
    end
    return (false, (ZZ(1), N))
end
################################################################################

function compute_ell_curve(N::fmpz, D::fmpz)
    R = ResidueRing(ZZ, N)
    if D == ZZ(-3)
        return (R.([0, -1]), ZZ(1))
    elseif D == ZZ(-4)
        return (R.([-1, 0]), ZZ(1))
    else
        T, t = ZZ["t"]
        p = hilbert_class_polynomial(Int(D), T)
        p = change_base_ring(R, p)
        j = roots(p)[1]
        c = inv(j - R(1728)) * j
        return ([-3c, -2c], ZZ(degree(p))) # Stores curve parameters and class number
    end
end

################################################################################

function random_point(v::Vector{fmpz_mod})
    R = parent(v[1])
    N = modulus(R)
    A = rand(Int, 1000)
    for a in A
        z = R(a)^3 + R(a) * v[1] + v[2]
        if jacobi_symbol(ZZ(z), N) != -1
            (bool, res) = issquare_with_square_root(z)
            if bool
                return (bool, [R(a), res, R(1)])
            end
        end
    end
    return (false, R.([0, 0, 0]))
end

################################################################################

function compute_p1p2(v::Vector{fmpz_mod}, m::fmpz, q::fmpz)
    (bool, P) = random_point(v)
    R = parent(v[1])
    !bool && return (false, (R.([0, 0, 0]), ZZ(1)), (R.([0, 0, 0]), ZZ(1)))
    P2 = _scalar_mult(P, v, div(m, q))
    P1 = _scalar_mult(P2[1], v, div(m, q))
    return (true, P1, P2)
end

################################################################################

function quadratic_non_residue(N::fmpz, D::fmpz)
    R = ResidueRing(ZZ, N)
    A = rand(Int, 1000)
    for a in A
        (bool, res) = issquare_with_square_root(R(a))
        if !bool
            if D == ZZ(-3)
                !isone(R(a)^div(N - 1, 3)) && return R(a)
            else
                return R(a)
            end
        end
    end
    return R(0)
end

################################################################################

function atkin_morain_step(
    N::fmpz,
    D::fmpz,
    m::fmpz,
    q::fmpz,
    res::Vector{Tuple{fmpz,Int}} = Vector{Tuple{fmpz,Int}}(),
    pos::Int = 1,
)
    (v, h) = compute_ell_curve(N, D)
    R = parent(v[1])
    O_E = R.([0, 1, 0])
    g = quadratic_non_residue(N, D)

    for k = 0:(h-1)
        P1 = (O_E, ZZ(1))
        while P1[1] == O_E
            (bool2, P1, P2) = compute_p1p2(v, m, q)
            !bool2 && return (false, N, ZZ(1), res, pos)
            P2[1] == R.([0, 0, 0]) && return (false, N, P2[2], res, pos)
            P1[1] == R.([0, 0, 0]) && return (false, N, P1[2], res, pos)
            P2[1] != O_E && return (true, q, ZZ(1), push!(res, (q, pos)), pos)
        end
        if D == ZZ(-3)
            v = [v[1], v[2] * g]
        elseif D == ZZ(-4)
            v = [v[1] * g, v[2]]
        else
            v = [v[1] * g^2, v[2] * g^3]
        end
    end
    return (false, N, ZZ(1), res, pos)
end

################################################################################

function atkin_morain(
    N::fmpz,
    Arr::Vector{fmpz} = Vector{fmpz}(),
    pos::Int = 1,
    res::Vector{Tuple{fmpz,Int}} = Vector{Tuple{fmpz,Int}}(),
)
    N < ZZ(10)^15 && return isprime(N)
    if length(Arr) == 0
        L = funddiscriminant(10^5)
    else
        L = Arr
    end
    pos > length(L) && error("No proof found")
    i = pos
    while i <= length(L)
        if jacobi_symbol(L[i], N) == 1
            (bool1, (m, q)) = find_q(N, L[i])
            if bool1
                (bool2, A, B, res, pos) = atkin_morain_step(N, L[i], m, q, res, i)
                if bool2
                    return atkin_morain(A, L, pos, res)
                elseif length(res) == 0
                    error("No proof found")
                else
                    # Backtracking
                    (A, pos) = res[length(res)]
                    return atkin_morain(A, L, pos, res)
                end
            end
        end
        i = i + 1
    end
    error("No proof found")
end

################################################################################
@doc Markdown.doc"""
    ECPP(n::fmpz)

The algorithm returns true if the number is prime, false if not, and an error if
it can't decide.
"""
function ECPP(n::fmpz)
    n == ZZ(2) && return true
    n % 2 == 0 && return false
    !Miller_Rabin_test(n) && return false
    L = funddiscriminant(10^4)
    return atkin_morain(n, L)
end
