export  f4, groebner_basis, groebner_basis_with_transformation_matrix, leading_ideal, syzygy_generators

# groebner stuff #######################################################
@doc Markdown.doc"""
    groebner_assure(I::MPolyIdeal, complete_reduction::Bool = false)
    groebner_assure(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool = false,
                    enforce_global_ordering::Bool = true)

**Note**: Internal function, subject to change, do not use.

Given an ideal `I` in a multivariate polynomial ring this function assures that a
Groebner basis w.r.t. the given monomial ordering is attached to `I` in `I.gb`.
It *currently* also ensures that the basis is defined on the Singular side in
`I.gb.S`, but this should not be relied upon: use `singular_assure(I.gb)` before
accessing `I.gb.S`.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ, ["x","y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*y-3*x,y^3-2*x^2*y])
ideal(x*y - 3*x, -2*x^2*y + y^3)

julia> Oscar.groebner_assure(I, degrevlex(gens(R)));

julia> I.gb[degrevlex(gens(R))]
Oscar.BiPolyArray{fmpq_mpoly}(Multivariate Polynomial Ring in x, y over Rational Field, #undef, Singular Polynomial Ring (QQ),(x,y),(dp(2),C), Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y - 3*x, y^3 - 6*x^2, 2*x^3 - 9*x), true, #undef, true)
```
"""
function groebner_assure(I::MPolyIdeal, complete_reduction::Bool = false)
    if isempty(I.gb)
        drl = default_ordering(base_ring(I))
        I.gb[drl] = groebner_assure(I, drl, complete_reduction)
        G = I.gb[drl]
    else
        G = first(values(I.gb))
    end
    return G
end

function groebner_assure(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool = false, enforce_global_ordering::Bool = true)
    return get!(I.gb, ordering) do
        _compute_groebner_basis(I.gens, ordering, complete_reduction, enforce_global_ordering)
    end
end

@doc Markdown.doc"""
    _compute_groebner_basis(B::BiPolyArray, ordering::MonomialOrdering,
                            complete_reduction::Bool = false, enforce_global_ordering::Bool = true)

**Note**: Internal function, subject to change, do not use.

Given an `BiPolyArray` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of the
ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`. The Groebner basis is then
returned in `B.S`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> A = Oscar.BiPolyArray([x*y-3*x,y^3-2*x^2*y])
Oscar.BiPolyArray{fmpq_mpoly}(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x*y - 3*x, -2*x^2*y + y^3], #undef, #undef, false, #undef, true)

julia> B = Oscar._compute_groebner_basis(A, degrevlex(gens(R)))
Oscar.BiPolyArray{fmpq_mpoly}(Multivariate Polynomial Ring in x, y over Rational Field, #undef, Singular Polynomial Ring (QQ),(x,y),(dp(2),C), Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y - 3*x, y^3 - 6*x^2, 2*x^3 - 9*x), true, #undef, true)
```
"""
function _compute_groebner_basis(B::BiPolyArray, ordering::MonomialOrdering, complete_reduction::Bool = false, enforce_global_ordering::Bool = true)
   singular_assure(B, ordering)
   R = B.Sx
   if enforce_global_ordering
     !Oscar.Singular.has_global_ordering(R) && error("The ordering has to be a global ordering.")
   end
   I  = Singular.Ideal(R, gens(B.S)...)
   i  = Singular.std(I, complete_reduction = complete_reduction)
   BA = BiPolyArray(B.Ox, i)
   BA.isGB  = true
   if isdefined(BA, :S)
       BA.S.isGB  = true
   end

   return BA
end

@doc Markdown.doc"""
    function groebner_basis(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool = false, enforce_global_ordering::Bool = true)

Given an ideal `I` and optional parameters monomial ordering `ordering` and `complete_reduction`,
compute a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of `I`
w.r.t. the given monomial ordering `ordering` (as default `degree reverse lexicographical`).
If `enforce_global_ordering = false` it is not checked whether the given ordering is global.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*y-3*x,y^3-2*x^2*y])
ideal(x*y - 3*x, -2*x^2*y + y^3)

julia> H = groebner_basis(I, ordering=lex(gens(R)))
3-element Vector{fmpq_mpoly}:
 y^4 - 3*y^3
 x*y - 3*x
 6*x^2 - y^3
```
"""
function groebner_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool=false, enforce_global_ordering::Bool = true)
    groebner_assure(I, ordering, complete_reduction, enforce_global_ordering)
    return collect(I.gb[ordering])
end

@doc Markdown.doc"""
    f4(I::MPolyIdeal, <keyword arguments>)

Compute a Groebner basis of the given ideal `I` w.r.t. to the degree reverse lexicographical monomial ordering using Faugère's F4 algorithm.
See [Fau99](@cite) for more information.

**Note**: At the moment only ground fields of characteristic `p`, `p` prime, `p < 2^{31}` are supported.

# Arguments
- `Ì::MPolyIdeal`: input ideal.
- `initial_hts::Int=17`: initial hash table size `log_2`.
- `nr_thrds::Int=1`: number of threads for parallel linear algebra.
- `max_nr_pairs::Int=0`: maximal number of pairs per matrix, only bounded by minimal degree if `0`.
- `la_option::Int=2`: linear algebra option: exact sparse-dense (`1`), exact sparse (`2`, default), probabilistic sparse-dense (`42`), probabilistic sparse(`44`).
- `eliminate::Int=0`: size of first block of variables to be eliminated.
- `complete_reduction::Bool=true`: compute a reduced Gröbner basis for `I`
- `info_level::Int=0`: info level printout: off (`0`, default), summary (`1`), detailed (`2`).

# Examples
```jldoctest
julia> R,(x,y,z) = PolynomialRing(GF(101), ["x","y","z"], ordering=:degrevlex)
(Multivariate Polynomial Ring in x, y, z over Galois field with characteristic 101, gfp_mpoly[x, y, z])

julia> I = ideal(R, [x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
ideal(x + 2*y + 2*z + 100, x^2 + 2*y^2 + 2*z^2 + 100*x, 2*x*y + 2*y*z + 100*y)

julia> f4(I)
4-element Vector{gfp_mpoly}:
 x + 2*y + 2*z + 100
 y*z + 82*z^2 + 10*y + 40*z
 y^2 + 60*z^2 + 20*y + 81*z
 z^3 + 28*z^2 + 64*y + 13*z

julia> isdefined(I, :gb)
true
```
"""
function f4(
        I::MPolyIdeal;
        initial_hts::Int=17,
        nr_thrds::Int=1,
        max_nr_pairs::Int=0,
        la_option::Int=2,
        eliminate::Int=0,
        complete_reduction::Bool=true,
        info_level::Int=0
        )
    AI = AlgebraicSolving.Ideal(I.gens.O)
    AlgebraicSolving.groebner_basis(AI,
                initial_hts = initial_hts,
                nr_thrds = nr_thrds,
                max_nr_pairs = max_nr_pairs,
                la_option = la_option,
                eliminate = eliminate,
                complete_reduction = complete_reduction,
                info_level = info_level)

    vars = gens(base_ring(I))[eliminate+1:end]
    I.gb[degrevlex(vars)] =
        BiPolyArray(AI.gb[eliminate], keep_ordering = false, isGB = true)

    return AI.gb[eliminate]
end

@doc Markdown.doc"""
    groebner_basis_with_transform(B::BiPolyArray, ordering::MonomialOrdering, complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `BiPolyArray` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of the
ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering` and the transformation matrix from the ideal to the Groebner basis. Return value is a BiPolyArray together with a map.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> A = Oscar.BiPolyArray([x*y-3*x,y^3-2*x^2*y])
Oscar.BiPolyArray{fmpq_mpoly}(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x*y - 3*x, -2*x^2*y + y^3], #undef, #undef, false, #undef, true)

julia> B,m = Oscar.groebner_basis_with_transform(A, degrevlex(gens(R)))
(Oscar.BiPolyArray{fmpq_mpoly}(Multivariate Polynomial Ring in x, y over Rational Field, #undef, Singular Polynomial Ring (QQ),(x,y),(dp(2),C), Singular ideal over Singular Polynomial Ring (QQ),(x,y),(dp(2),C) with generators (x*y - 3*x, y^3 - 6*x^2, 6*x^3 - 27*x), false, #undef, true), [1 2*x -2*x^2+y^2+3*y+9; 0 1 -x])
```
"""
function groebner_basis_with_transform(B::BiPolyArray, ordering::MonomialOrdering, complete_reduction::Bool = false)
   if !isdefined(B, :ordering)
      singular_assure(B, ordering)
   elseif ordering != B.ordering
     R = singular_poly_ring(B.Ox, ordering)
     i = Singular.Ideal(R, [R(x) for x = B])
     i, m = Singular.lift_std(i, complete_reduction = complete_reduction)
     return BiPolyArray(B.Ox, i), map_entries(x->B.Ox(x), m)
   end

   if !isdefined(B, :S)
     B.S = Singular.Ideal(B.Sx, [B.Sx(x) for x = B.O])
   end

   i, m = Singular.lift_std(B.S, complete_reduction = complete_reduction)
   return BiPolyArray(B.Ox, i), map_entries(x->B.Ox(x), m)
 end

@doc Markdown.doc"""
    groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool=false)

Return a pair `G, m` where `G` is a Groebner basis of the ideal `I` with respect to the
monomial ordering `ordering` (as default degree reverse lexicographical), and `m` is a transformation matrix from
gens(I)` to `G`. If complete_reduction`is set to `true` (as default 'false') then `G` will be the reduced Groebner basis.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ,["x","y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*y^2-1,x^3+y^2+x*y])
ideal(x*y^2 - 1, x^3 + x*y + y^2)

julia> G,m = groebner_basis_with_transformation_matrix(I)
(fmpq_mpoly[x*y^2 - 1, x^3 + x*y + y^2, x^2 + y^4 + y], [1 0 -x^2-y; 0 1 y^2])

julia> gens(I)*m == G
true

```
"""
function groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool = false)
   G, m = groebner_basis_with_transform(I.gens, ordering, complete_reduction)
   I.gb[ordering]  = G
   return collect(G), m
 end

# syzygies #######################################################
@doc Markdown.doc"""
    syzygy_generators(a::Vector{<:MPolyElem})

Return generators for the syzygies on the given polynomials.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> S = syzygy_generators([x^3+y+2,x*y^2-13*x^2,y-14])
3-element Vector{FreeModElem{fmpq_mpoly}}:
 (-y + 14)*e[2] + (-13*x^2 + x*y^2)*e[3]
 (-169*y + 2366)*e[1] + (-13*x*y + 182*x - 196*y + 2744)*e[2] + (13*x^2*y^2 - 2548*x^2 + 196*x*y^2 + 169*y + 338)*e[3]
 (-13*x^2 + 196*x)*e[1] + (-x^3 - 16)*e[2] + (x^4*y + 14*x^4 + 13*x^2 + 16*x*y + 28*x)*e[3]
```
"""
function syzygy_generators(a::Vector{<:MPolyElem})
  I = ideal(a)
  singular_assure(I)
  s = Singular.syz(I.gens.S)
  F = free_module(parent(a[1]), length(a))
  @assert rank(s) == length(a)
  return [F(s[i]) for i=1:Singular.ngens(s)]
end

# leading ideal #######################################################
@doc Markdown.doc"""
    leading_ideal(g::Vector{T}; ordering::MonomialOrdering) where { T <: MPolyElem }

Return the ideal generated by the leading monomials of the given polynomials
w.r.t. the given monomial ordering.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5], ordering=degrevlex(gens(R)))
ideal(x*y^2, y^5)

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5], ordering=lex(gens(R)))
ideal(x*y^2, x^3)
```
"""
function leading_ideal(g::Vector{T}; ordering::MonomialOrdering) where { T <: MPolyElem }
    return ideal(parent(g[1]), [first(monomials(f, ordering)) for f in g])
end


@doc Markdown.doc"""
    leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering)

Given a multivariate polynomial ideal `Ì` this function returns the
leading ideal for `I`. This is done w.r.t. the given monomial ordering.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R,[x*y^2-3*x, x^3-14*y^5])
ideal(x*y^2 - 3*x, x^3 - 14*y^5)

julia> L = leading_ideal(I, ordering=degrevlex(gens(R)))
ideal(x*y^2, x^4, y^5)

julia> L = leading_ideal(I, ordering=lex(gens(R)))
ideal(y^7, x*y^2, x^3)
```
"""
function leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering)
  G = groebner_basis(I, ordering=ordering, enforce_global_ordering=false)
  return ideal(base_ring(I), [first(monomials(g, ordering)) for g in G])
end

@doc Markdown.doc"""
    normal_form_internal(I::Singular.sideal, J::MPolyIdeal)

**Note**: Internal function, subject to change, do not use.

Compute the normal form of the generators `gens(I)` of the ideal `I` w.r.t. a
Groebner basis of `J`.

CAVEAT: This computation needs a Groebner basis of `J`. If this Groebner basis
is not available, one is computed automatically. This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> groebner_basis(J)
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> SR = singular_poly_ring(base_ring(J))
Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C)

julia> I = Singular.Ideal(SR,[SR(-1+c+b+a^3),SR(-1+b+c*a+2*a^3),SR(5+c*b+c^2*a)])
Singular ideal over Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C) with generators (a^3 + b + c - 1, 2*a^3 + a*c + b - 1, a*c^2 + b*c + 5)

julia> Oscar.normal_form_internal(I,J)
3-element Vector{fmpq_mpoly}:
 a^3
 2*a^3 + 2*a - 2*c
 4*a - 2*c^2 - c + 5
```
"""
function normal_form_internal(I::Singular.sideal, J::MPolyIdeal)
  groebner_assure(J)
  G = collect(values(J.gb))[1]
  singular_assure(G)
  K = ideal(base_ring(J), reduce(I, G.S))
  return [J.gens.Ox(x) for x = gens(K.gens.S)]
end

@doc Markdown.doc"""
    normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }

Compute the normal form of the polynomial `f` w.r.t. a
Groebner basis of `J`.

CAVEAT: This computation needs a Groebner basis of `J`. If this Groebner basis
is not available, one is computed automatically. This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> groebner_basis(J)
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(-1+c+b+a^3, J)
a^3
```
"""
function normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }
    singular_assure(J)
    I = Singular.Ideal(J.gens.Sx, J.gens.Sx(f))
    N = normal_form_internal(I, J)
    return N[1]
end

@doc Markdown.doc"""
    normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }

Compute the normal form of the elements of the array `A` w.r.t. a
Groebner basis of `J`.

CAVEAT: This computation needs a Groebner basis of `J`. If this Groebner basis
is not available, one is computed automatically. This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> A = [-1+c+b+a^3,-1+b+c*a+2*a^3,5+c*b+c^2*a]
3-element Vector{fmpq_mpoly}:
 a^3 + b + c - 1
 2*a^3 + a*c + b - 1
 a*c^2 + b*c + 5

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> groebner_basis(J)
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(A, J)
3-element Vector{fmpq_mpoly}:
 a^3
 2*a^3 + 2*a - 2*c
 4*a - 2*c^2 - c + 5
```
"""
function normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }
    singular_assure(J)
    I = Singular.Ideal(J.gens.Sx, [J.gens.Sx(x) for x in A])
    normal_form_internal(I, J)
end

# standard basis for non-global orderings #############################
@doc Markdown.doc"""
    std_basis(I::MPolyIdeal, o::MonomialOrdering)

Compute a standard basis of `I` for the monomial ordering `o`.

**Note:** Since there is no general notion of complete reduction for 
non-global orderings, this option is not available for this command;
instead, use `groebner_basis` directly. 

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ, ["x","y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal([x*(x+1), x^2-y^2+(x-2)*y])
ideal(x^2 + x, x^2 + x*y - y^2 - 2*y)

julia> std_basis(I, negdegrevlex(gens(R)))
2-element Vector{fmpq_mpoly}:
 x
 y
```
"""
function std_basis(I::MPolyIdeal, o::MonomialOrdering)
  return groebner_basis(I, ordering=o, enforce_global_ordering=false, 
                        complete_reduction=false)
end

function normal_form(f::MPolyElem, J::MPolyIdeal, o::MonomialOrdering)
  groebner_assure(J, o, false, false)
  stdJ = J.gb[o]
  Sx = stdJ.Sx
  Ox = parent(f)
  I = Singular.Ideal(Sx, Sx(f))
  return Ox(gens(Singular.reduce(I, stdJ.S))[1])
end
