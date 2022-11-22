export  reduce, reduce_with_quotients, reduce_with_quotients_and_unit, f4,
		standard_basis, groebner_basis, standard_basis_with_transformation_matrix,
		groebner_basis_with_transformation_matrix,
		leading_ideal, syzygy_generators, is_standard_basis, is_groebner_basis

# groebner stuff #######################################################
@doc Markdown.doc"""
    groebner_assure(I::MPolyIdeal, complete_reduction::Bool = false)
    groebner_assure(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool = false)

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

julia> Oscar.groebner_assure(I, degrevlex(R));

julia> I.gb[degrevlex(R)]
Gröbner basis with elements
1 -> x*y - 3*x
2 -> y^3 - 6*x^2
3 -> 2*x^3 - 9*x
with respect to the ordering
degrevlex([x, y])
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

function groebner_assure(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool = false)
    return get!(I.gb, ordering) do
        _compute_standard_basis(I.gens, ordering, complete_reduction)
    end
end

@doc Markdown.doc"""
    _compute_standard_basis(B::IdealGens; ordering::MonomialOrdering,
                            complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `IdealGens` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a Groebner basis (if `complete_reduction = true` the reduced Groebner basis) of the
ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`. The Groebner basis is then
returned in `B.S`.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ, ["x","y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> A = Oscar.IdealGens([x*y-3*x,y^3-2*x^2*y])
Ideal generating system with elements
1 -> x*y - 3*x
2 -> -2*x^2*y + y^3

julia> B = Oscar._compute_standard_basis(A, degrevlex(R))
Gröbner basis with elements
1 -> x*y - 3*x
2 -> y^3 - 6*x^2
3 -> 2*x^3 - 9*x
with respect to the ordering
degrevlex([x, y])
```
"""
function _compute_standard_basis(B::IdealGens, ordering::MonomialOrdering, complete_reduction::Bool = false)
	singular_assure(B, ordering)
	R = B.Sx
	I  = Singular.Ideal(R, gens(B.S)...)
	i  = Singular.std(I, complete_reduction = complete_reduction)
	BA = IdealGens(B.Ox, i, complete_reduction)
	BA.isGB = true
	BA.ord = ordering
	if isdefined(BA, :S)
	   BA.S.isGB  = true
	end
	return BA
end

# standard basis for non-global orderings #############################
@doc Markdown.doc"""
    standard_basis(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool = false)

Return a standard basis of `I` with respect to `ordering`. 

!!! note
    The returned standard basis is reduced if `ordering` is `global` and `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ, ["x","y"]);

julia> I = ideal([x*(x+1), x^2-y^2+(x-2)*y]);

julia> standard_basis(I, ordering = negdegrevlex(R))
Standard basis with elements
1 -> x
2 -> y
with respect to the ordering
negdegrevlex([x, y])
```
"""
function standard_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)),
	complete_reduction::Bool = false, algorithm::Symbol = :direct)
	complete_reduction && @assert is_global(ordering)
	if algorithm == :direct
		if !haskey(I.gb, ordering)
			I.gb[ordering] = _compute_standard_basis(I.gens, ordering, complete_reduction)
		end
		if !I.gb[ordering].isReduced
			I.gb[ordering] = _compute_standard_basis(I.gb[ordering].gens, ordering, complete_reduction)
		end
	elseif algorithm == :fglm
		_compute_groebner_basis_using_fglm(I, ordering, default_ordering(base_ring(I)))
	elseif algorithm == :f4
		_compute_groebner_basis_using_f4(I, complete_reduction=complete_reduction)
	end
	return I.gb[ordering]
end

@doc Markdown.doc"""
    groebner_basis(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool = false)

If `ordering` is global, return a Gröbner basis of `I` with respect to `ordering`.

!!! note
    The returned Gröbner basis is reduced if `complete_reduction = true`.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> G = groebner_basis(I)
Gröbner basis with elements
1 -> y^2 - x*z
2 -> x*y - z
3 -> x^2 - y
with respect to the ordering
degrevlex([x, y, z])

julia> elements(G)
3-element Vector{fmpq_mpoly}:
 -x*z + y^2
 x*y - z
 x^2 - y

julia> elements(G) == gens(G)
true

julia> groebner_basis(I, ordering = lex(R))
Gröbner basis with elements
1 -> y^3 - z^2
2 -> x*z - y^2
3 -> x*y - z
4 -> x^2 - y
with respect to the ordering
lex([x, y, z])
```
```jldoctest
julia> R, (x, y) = GradedPolynomialRing(QQ, ["x", "y"], [1, 3]);

julia> I = ideal(R, [x*y-3*x^4,y^3-2*x^6*y]);

julia> groebner_basis(I)
Gröbner basis with elements
1 -> 3*x^4 - x*y
2 -> 2*x^3*y^2 - 3*y^3
3 -> x*y^3
4 -> y^4
with respect to the ordering
wdegrevlex([x, y], [1, 3])
```
"""
function groebner_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool=false,
	algorithm::Symbol = :direct)
    is_global(ordering) || error("Ordering must be global")
    return standard_basis(I, ordering=ordering, complete_reduction=complete_reduction, algorithm=algorithm)
end

@doc Markdown.doc"""
    f4(I::MPolyIdeal, <keyword arguments>)

Compute a Gröbner basis of `I` with respect to `degrevlex` using Faugère's F4 algorithm.
See [Fau99](@cite) for more information.

!!! note
    At current state only prime fields of characteristic `0 < p < 2^{31}` are supported.

# Possible keyword arguments
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
Gröbner basis with elements
1 -> x + 2*y + 2*z + 100
2 -> y*z + 82*z^2 + 10*y + 40*z
3 -> y^2 + 60*z^2 + 20*y + 81*z
4 -> z^3 + 28*z^2 + 64*y + 13*z
with respect to the ordering
degrevlex([x, y, z])

julia> isdefined(I, :gb)
true
```
"""
function _compute_groebner_basis_using_f4(
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
    ord = degrevlex(vars)
    I.gb[ord] =
        IdealGens(AI.gb[eliminate], ord, keep_ordering = false, isGB = true)
	I.gb[ord].isReduced = complete_reduction

    return I.gb[ord]
end

@doc Markdown.doc"""
    _compute_standard_basis_with_transform(B::BiPolyArray, ordering::MonomialOrdering, complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `IdealGens` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a standard basis (if `ordering` is a global monomial ordering and `complete_reduction = true`
the reduced Groebner basis) of the ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`
and the transformation matrix from the ideal to the standard basis. Return value is a IdealGens together with a map.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> A = Oscar.IdealGens([x*y-3*x,y^3-2*x^2*y])
Ideal generating system with elements
1 -> x*y - 3*x
2 -> -2*x^2*y + y^3

julia> B,m = Oscar._compute_standard_basis_with_transform(A, degrevlex(R))
(Ideal generating system with elements
1 -> x*y - 3*x
2 -> -6*x^2 + y^3
3 -> 6*x^3 - 27*x, [1 2*x -2*x^2+y^2+3*y+9; 0 1 -x])
```
"""
function _compute_standard_basis_with_transform(B::IdealGens, ordering::MonomialOrdering, complete_reduction::Bool = false)
   if !isdefined(B, :ordering)
      singular_assure(B, ordering)
   elseif ordering != B.ordering
     R = singular_poly_ring(B.Ox, ordering)
     i = Singular.Ideal(R, [R(x) for x = B])
     i, m = Singular.lift_std(i, complete_reduction = complete_reduction)
     return IdealGens(B.Ox, i), map_entries(x->B.Ox(x), m)
   end

   if !isdefined(B, :S)
     B.S = Singular.Ideal(B.Sx, [B.Sx(x) for x = B.O])
   end

   i, m = Singular.lift_std(B.S, complete_reduction = complete_reduction)
   return IdealGens(B.Ox, i), map_entries(x->B.Ox(x), m)
 end

@doc Markdown.doc"""
    standard_basis_with_transformation_matrix(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool=false)

Return a pair `G`, `T`, say, where `G` is a standard basis of `I` with respect to `ordering`, and `T` 
is a transformation matrix from `gens(I)` to `G`. That is, `gens(I)*T == G`.

!!! note
    The returned Gröbner basis is reduced if `ordering` is a global monomial odering and `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ,["x","y"]);

julia> I = ideal([x*y^2-1,x^3+y^2+x*y]);

julia> G, T = standard_basis_with_transformation_matrix(I, ordering=neglex(R))
(Standard basis with elements
1 -> 1 - x*y^2
with respect to the ordering
neglex([x, y]), [-1; 0])

julia> gens(I)*T == gens(G)
true
```
"""
function standard_basis_with_transformation_matrix(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool = false)
	complete_reduction && @assert is_global(ordering)
	G, m = _compute_standard_basis_with_transform(I.gens, ordering, complete_reduction)
	G.isGB = true
	I.gb[ordering]  = G
	return G, m
 end

@doc Markdown.doc"""
    groebner_basis_with_transformation_matrix(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool=false)

Return a pair `G`, `T`, say, where `G` is a Gröbner basis of `I` with respect to `ordering`, and `T` 
is a transformation matrix from `gens(I)` to `G`. That is, `gens(I)*T == G`.

!!! note
    The returned Gröbner basis is reduced if `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = PolynomialRing(QQ,["x","y"]);

julia> I = ideal([x*y^2-1,x^3+y^2+x*y]);

julia> G, T = groebner_basis_with_transformation_matrix(I)
(Gröbner basis with elements
1 -> x*y^2 - 1
2 -> x^3 + x*y + y^2
3 -> y^4 + x^2 + y
with respect to the ordering
degrevlex([x, y]), [1 0 -x^2-y; 0 1 y^2])

julia> gens(I)*T == gens(G)
true
```
"""
function groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool = false)
    is_global(ordering) || error("Ordering must be global")
	return standard_basis_with_transformation_matrix(I, ordering=ordering, complete_reduction=complete_reduction)
 end

# syzygies #######################################################
@doc Markdown.doc"""
    syzygy_generators(G::Vector{<:MPolyElem})

Return generators for the syzygies on the polynomials given as elements of `G`.

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
leading_ideal(G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(G[1]))) where { T <: MPolyElem }

Return the leading ideal of `G` with respect to `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5], ordering=degrevlex(R))
ideal(x*y^2, y^5)

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5], ordering=lex(R))
ideal(x*y^2, x^3)
```
"""
function leading_ideal(G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(G[1]))) where { T <: MPolyElem }
    return ideal(parent(G[1]), [leading_monomial(f; ordering = ordering) for f in G])
end

function leading_ideal(I::IdealGens{T}) where { T <: MPolyElem }
    return ideal(base_ring(I), [leading_monomial(f; ordering = I.ord) for f in I])
end

function leading_ideal(I::IdealGens{T}, ordering::MonomialOrdering) where {T <: MPolyElem}
    return ideal(base_ring(I), [leading_monomial(f; ordering = ordering) for f in I])
end


@doc Markdown.doc"""
    leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)))

Return the leading ideal of `I` with respect to `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"])
(Multivariate Polynomial Ring in x, y over Rational Field, fmpq_mpoly[x, y])

julia> I = ideal(R,[x*y^2-3*x, x^3-14*y^5])
ideal(x*y^2 - 3*x, x^3 - 14*y^5)

julia> L = leading_ideal(I, ordering=degrevlex(R))
ideal(x*y^2, x^4, y^5)

julia> L = leading_ideal(I, ordering=lex(R))
ideal(y^7, x*y^2, x^3)
```
"""
function leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)))
  G = groebner_basis(I, ordering=ordering)
  return ideal(base_ring(I), [leading_monomial(g; ordering = ordering) for g in G])
end

@doc Markdown.doc"""
    normal_form_internal(I::Singular.sideal, J::MPolyIdeal, o::MonomialOrdering)

**Note**: Internal function, subject to change, do not use.

Compute the normal form of the generators `gens(I)` of the ideal `I` w.r.t. a
Groebner basis of `J` and the monomial ordering `o`.

CAVEAT: This computation needs a Groebner basis of `J` and the monomial ordering
`o`. If this Groebner basis is not available, one is computed automatically.
This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> gens(groebner_basis(J))
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> SR = singular_poly_ring(base_ring(J))
Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C)

julia> I = Singular.Ideal(SR,[SR(-1+c+b+a^3),SR(-1+b+c*a+2*a^3),SR(5+c*b+c^2*a)])
Singular ideal over Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C) with generators (a^3 + b + c - 1, 2*a^3 + a*c + b - 1, a*c^2 + b*c + 5)

julia> Oscar.normal_form_internal(I,J,default_ordering(base_ring(J)))
3-element Vector{fmpq_mpoly}:
 a^3
 2*a^3 + 2*a - 2*c
 4*a - 2*c^2 - c + 5
```
"""
function normal_form_internal(I::Singular.sideal, J::MPolyIdeal, o::MonomialOrdering)
  groebner_assure(J, o)
  G = J.gb[o]  
  singular_assure(G, o)
  K = ideal(base_ring(J), reduce(I, G.S))
  return [J.gens.Ox(x) for x = gens(K.gens.S)]
end

@doc Markdown.doc"""
	reduce(I::IdealGens, J::IdealGens; 
          ordering::MonomialOrdering = default_ordering(base_ring(J)))

Return a `Vector` whose elements are the underlying elements of `I`
reduced by the underlying generators of `J` w.r.t. the monomial
ordering `ordering`. `J` need not be a Groebner basis. The returned
`Vector` will have the same number of elements as `I`, even if they
are zero.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(GF(11), ["x", "y", "z"]);

julia> I = ideal(R, [x^2, x*y - y^2]);

julia> J = ideal(R, [y^3])
ideal(y^3)

julia> reduce(J.gens, I.gens)
1-element Vector{gfp_mpoly}:
 y^3

julia> reduce(J.gens, groebner_basis(I))
1-element Vector{gfp_mpoly}:
 0

julia> reduce(y^3, [x^2, x*y-y^3])
x*y

julia> reduce(y^3, [x^2, x*y-y^3], ordering=lex(R))
y^3

julia> reduce([y^3], [x^2, x*y-y^3], ordering=lex(R))
1-element Vector{gfp_mpoly}:
 y^3
```
"""
function reduce(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))
	@assert base_ring(J) == base_ring(I)
	singular_assure(I, ordering)
	singular_assure(J, ordering)
	res = reduce(I.gens.S, J.gens.S)
	return [J.gens.Ox(x) for x = gens(res)]
end

@doc Markdown.doc"""
	reduce(g::T, F::Vector{T}; 
          ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}

If `ordering` is global, return the remainder in a standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.
Otherwise, return the remainder in a *weak* standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.

	reduce(G::Vector{T}, F::Vector{T};
          ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}

Return a `Vector` which contains, for each element `g` of `G`, a remainder as above.

!!! note
    In the global case, the returned remainders are fully reduced.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(QQ, ["x", "y"]);

julia> reduce(y^3, [x^2, x*y-y^3])
x*y

julia> reduce(y^3, [x^2, x*y-y^3], ordering = lex(R))
y^3
```

```jldoctest
julia> R, (z, y, x) = PolynomialRing(QQ, ["z", "y", "x"]);

julia> f1 = y-x^2; f2 = z-x^3;

julia> g = x^3*y-3*y^2*z^2+x*y*z;

julia> reduce(g, [f1, f2], ordering = lex(R))
-3*x^10 + x^6 + x^5
```
"""
function reduce(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(f))) where {T <: MPolyElem}
	@assert parent(f) == parent(F[1])
	R = parent(f)
	I = IdealGens(R, [f], ordering)
	J = IdealGens(R, F, ordering)
	redv = reduce(I, J, ordering=ordering)
	return redv[1]
end

function reduce(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
	@assert parent(F[1]) == parent(G[1])
	R = parent(F[1])
	I = IdealGens(R, F, ordering)
	J = IdealGens(R, G, ordering)
	return reduce(I, J, ordering=ordering)
end

@doc Markdown.doc"""
	reduce_with_quotients_and_unit(g::T, F::Vector{T};
          ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}

Return the unit, the quotients and the remainder in a weak standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.

	reduce_with_quotients(G::Vector{T}, F::Vector{T};
          ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}

Return a `Vector` which contains, for each element `g` of `G`, a unit, quotients, and a remainder as above.

!!! note
    In the global case, a standard representation with a fully reduced remainder is computed.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);

julia> f1 = x^2+x^2*y; f2 = y^3+x*y*z; f3 = x^3*y^2+z^4;

julia> g = x^3*y+x^5+x^2*y^2*z^2+z^6;

julia> u, Q, h = reduce_with_quotients_and_unit(g, [f1,f2, f3], ordering = negdegrevlex(R))
([y+1], [x^3-x*y^2*z^2+x*y+y^2*z^2 0 y*z^2+z^2], 0)

julia> u*g == Q[1]*f1+Q[2]*f2+Q[3]*f3+h
true

julia> G = [g, x*y^3-3*x^2*y^2*z^2];

julia> U, Q,  H = reduce_with_quotients_and_unit(G, [f1, f2, f3], ordering = lex(R));

julia> U
[1   0]
[0   1]

julia> H
2-element Vector{fmpq_mpoly}:
 -z^9 + z^7 + z^6 + z^4
 -3*z^7 + z^6

julia> U*G == Q*[f1, f2, f3]+H
true
```
"""
function reduce_with_quotients_and_unit(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
	@assert parent(f) == parent(F[1])
	R = parent(f)
	I = IdealGens(R, [f], ordering)
	J = IdealGens(R, F, ordering)
	u, q, r = _reduce_with_quotients_and_unit(I, J, ordering)
	return u, q, r[1]
end

function reduce_with_quotients_and_unit(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
	@assert parent(F[1]) == parent(G[1])
	R = parent(F[1])
	I = IdealGens(R, F, ordering)
	J = IdealGens(R, G, ordering)
	return _reduce_with_quotients_and_unit(I, J, ordering)
end

@doc Markdown.doc"""
        reduce_with_quotients_and_unit(I::IdealGens, J::IdealGens; 
          ordering::MonomialOrdering = default_ordering(base_ring(J)))

Return a `Tuple` consisting of a `Generic.MatSpaceElem` `M`, a
`Vector` `res` whose elements are the underlying elements of `I`
reduced by the underlying generators of `J` w.r.t. the monomial
ordering `ordering` and a diagonal matrix `units` such that `M *
gens(J) + res == units * gens(I)`. If `ordering` is global then
`units` will always be the identity matrix, see also
`reduce_with_quotients`. `J` need not be a Groebner basis. `res` will
have the same number of elements as `I`, even if they are zero.

# Examples
```jldoctest
julia> R, (x, y) = PolynomialRing(GF(11), ["x", "y"]);

julia> I = ideal(R, [x]);

julia> R, (x, y) = PolynomialRing(GF(11), ["x", "y"]);

julia> I = ideal(R, [x]);

julia> J = ideal(R, [x+1]);

julia> unit, M, res = reduce_with_quotients_and_unit(I.gens, J.gens, ordering = neglex(R))
([x+1], [x], gfp_mpoly[0])

julia> M * gens(J) + res == unit * gens(I)
true

julia> f = x^3*y^2-y^4-10
x^3*y^2 + 10*y^4 + 1

julia> F = [x^2*y-y^3, x^3-y^4]
2-element Vector{gfp_mpoly}:
 x^2*y + 10*y^3
 x^3 + 10*y^4

julia> reduce_with_quotients_and_unit(f, F)
([1], [x*y 10*x+1], x^4 + 10*x^3 + 1)

julia> unit, M, res = reduce_with_quotients_and_unit(f, F, ordering=lex(R))
([1], [0 y^2], y^6 + 10*y^4 + 1)

julia> M * F + [res] == unit * [f]
true
```
"""
function reduce_with_quotients_and_unit(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))
	return _reduce_with_quotients_and_unit(I, J, ordering)
end


@doc Markdown.doc"""
        reduce_with_quotients(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))

Return a `Tuple` consisting of a `Generic.MatSpaceElem` `M` and a
`Vector` `res` whose elements are the underlying elements of `I`
reduced by the underlying generators of `J` w.r.t. the monomial
ordering `ordering` such that `M * gens(J) + res == gens(I)` if `ordering` is global.
If `ordering` is local then this equality holds after `gens(I)` has been multiplied
with an unkown diagonal matrix of units, see reduce_with_quotients_and_unit` to
obtain this matrix. `J` need not be a Groebner basis. `res` will have the same number
of elements as `I`, even if they are zero.

# Examples
```jldoctest
julia> R, (x, y, z) = PolynomialRing(GF(11), ["x", "y", "z"]);

julia> J = ideal(R, [x^2, x*y - y^2]);

julia> I = ideal(R, [x*y, y^3]);

julia> gb = groebner_basis(J)
Gröbner basis with elements
1 -> x*y + 10*y^2
2 -> x^2
3 -> y^3
with respect to the ordering
degrevlex([x, y, z])

julia> M, res = reduce_with_quotients(I.gens, gb)
([1 0 0; 0 0 1], gfp_mpoly[y^2, 0])

julia> M * gens(gb) + res == gens(I)
true

julia> f = x^3*y^2-y^4-10
x^3*y^2 + 10*y^4 + 1

julia> F = [x^2*y-y^3, x^3-y^4]
2-element Vector{gfp_mpoly}:
 x^2*y + 10*y^3
 x^3 + 10*y^4

julia> reduce_with_quotients_and_unit(f, F)
([1], [x*y 10*x+1], x^4 + 10*x^3 + 1)

julia> unit, M, res = reduce_with_quotients_and_unit(f, F, ordering=lex(R))
([1], [0 y^2], y^6 + 10*y^4 + 1)

julia> M * F + [res] == unit * [f]
true
```
"""
function reduce_with_quotients(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))
    _, q, r = _reduce_with_quotients_and_unit(I, J, ordering)
    return q, r
end

@doc Markdown.doc"""
	reduce_with_quotients(g::T, F::Vector{T}; 
          ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}

If `ordering` is global, return the quotients and the remainder in a standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.
Otherwise, return the quotients and the remainder in a *weak* standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.

	reduce_with_quotients(G::Vector{T}, F::Vector{T}; 
          ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}

Return a `Vector` which contains, for each element `g` of `G`, quotients and a remainder as above.

!!! note
    In the global case, the returned remainders are fully reduced.

# Examples

```jldoctest
julia> R, (z, y, x) = PolynomialRing(QQ, ["z", "y", "x"]);

julia> f1 = y-x^2; f2 = z-x^3;

julia> g = x^3*y-3*y^2*z^2+x*y*z;

julia> Q, h = reduce_with_quotients(g, [f1, f2], ordering = lex(R));

julia> Q
[-3*y*x^6 - 3*x^8 + x^4 + x^3   -3*z*y^2 - 3*y^2*x^3 + y*x]

julia> h
-3*x^10 + x^6 + x^5

julia> g == Q[1]*f1+Q[2]*f2+h
true

julia> G = [g, x*y^3-3*x^2*y^2*z^2];

julia> Q, H = reduce_with_quotients(G, [f1, f2], ordering = lex(R));

julia> Q
[          -3*y*x^6 - 3*x^8 + x^4 + x^3   -3*z*y^2 - 3*y^2*x^3 + y*x]
[y^2*x - 3*y*x^8 + y*x^3 - 3*x^10 + x^5     -3*z*y^2*x^2 - 3*y^2*x^5]

julia> H
2-element Vector{fmpq_mpoly}:
 -3*x^10 + x^6 + x^5
 -3*x^12 + x^7

julia> G == Q*[f1, f2]+H
true
```
"""
function reduce_with_quotients(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
	@assert parent(f) == parent(F[1])
	R = parent(f)
	I = IdealGens(R, [f], ordering)
	J = IdealGens(R, F, ordering)
	_, q, r = _reduce_with_quotients_and_unit(I, J, ordering)
	return q, r[1]
end

function reduce_with_quotients(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyElem}
	@assert parent(F[1]) == parent(G[1])
	R = parent(F[1])
	I = IdealGens(R, F, ordering)
	J = IdealGens(R, G, ordering)
	_, q, r = _reduce_with_quotients_and_unit(I, J, ordering)
	return q, r
end

function _reduce_with_quotients_and_unit(I::IdealGens, J::IdealGens, ordering::MonomialOrdering = default_ordering(base_ring(J)))
	@assert base_ring(J) == base_ring(I)
	singular_assure(I, ordering)
	singular_assure(J, ordering)
	res = Singular.division(I.gens.S, J.gens.S)
	return matrix(base_ring(I), res[3]), matrix(base_ring(I), res[1]), [J.gens.Ox(x) for x = gens(res[2])]
end

@doc Markdown.doc"""
    normal_form(g::T, I::MPolyIdeal; 
      ordering::MonomialOrdering = default_ordering(base_ring(I))) where { T <: MPolyElem }

Compute the normal form of `g` mod `I` with respect to `ordering`.

    normal_form(G::Vector{T}, I::MPolyIdeal; 
      ordering::MonomialOrdering=default_ordering(base_ring(I))) where { T <: MPolyElem }

Return a `Vector` which contains for each element `g` of `G` a normal form as above.

# Examples
```jldoctest
julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> gens(groebner_basis(J))
2-element Vector{fmpq_mpoly}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(-1+c+b+a^3, J)
a^3

julia> R,(a,b,c) = PolynomialRing(QQ,["a","b","c"])
(Multivariate Polynomial Ring in a, b, c over Rational Field, fmpq_mpoly[a, b, c])

julia> A = [-1+c+b+a^3,-1+b+c*a+2*a^3,5+c*b+c^2*a]
3-element Vector{fmpq_mpoly}:
 a^3 + b + c - 1
 2*a^3 + a*c + b - 1
 a*c^2 + b*c + 5

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> gens(groebner_basis(J))
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
function normal_form(f::T, J::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(J))) where { T <: MPolyElem }
    singular_assure(J, ordering)
    I = Singular.Ideal(J.gens.Sx, J.gens.Sx(f))
    N = normal_form_internal(I, J, ordering)
    return N[1]
end

function normal_form(A::Vector{T}, J::MPolyIdeal; ordering::MonomialOrdering=default_ordering(base_ring(J))) where { T <: MPolyElem }
    singular_assure(J, ordering)
    I = Singular.Ideal(J.gens.Sx, [J.gens.Sx(x) for x in A])
    normal_form_internal(I, J, ordering)
end


function is_standard_basis(F::IdealGens; ordering::MonomialOrdering=default_ordering(base_ring(F)))
	if F.isGB && F.ord == ordering
		return true
	else
		# Try to reduce all possible s-polynomials, i.e. Buchberger's criterion
		R = base_ring(F)
		for i in 1:length(F)
			lt_i = leading_term(F[i], ordering=ordering)
			for j in i+1:length(F)
				lt_j = leading_term(F[j], ordering=ordering)
				lcm_ij  = lcm(lt_i, lt_j)
				sp_ij = div(lcm_ij, lt_i) * F[i] - div(lcm_ij, lt_j) * F[j]
				if reduce(IdealGens([sp_ij], ordering), F, ordering=ordering) != [R(0)]
					return false
				end
			end
		end
		F.isGB = true
		F.ord = ordering
		return true
	end
end

function is_groebner_basis(F::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(F)))
    is_global(ordering) || error("Ordering must be global")
	return is_standard_basis(F, ordering=ordering)
end

@doc Markdown.doc"""


Converts a Groebner basis `G` w.r.t. a given monomial ordering `<_1` for some ideal `I`
to a Groebner basis `H` w.r.t. another monomial ordering `<_2` for `I`.
"""
function fglm(G::IdealGens; ordering::MonomialOrdering)
	@assert G.isGB == true 
	@assert G.isReduced == true
	singular_assure(G)
	SR_destination, = Singular.PolynomialRing(base_ring(G.Sx),["$i" for i in gens(G.Sx)]; ordering = Singular.ordering_as_symbol(singular(ordering)))

	ptr = Singular.libSingular.fglmzero(G.S.ptr, G.Sx.ptr, SR_destination.ptr)
	return IdealGens(base_ring(G), Singular.sideal{Singular.spoly}(SR_destination, ptr, true))
end

function _compute_groebner_basis_using_fglm(I::MPolyIdeal,
	destination_ordering::MonomialOrdering, start_ordering::MonomialOrdering = default_ordering(base_ring(I)))
	haskey(I.gb, destination_ordering) && return I.gb[destination_ordering]
	@assert is_global(start_ordering) && is_global(destination_ordering)
	if start_ordering == degrevlex(base_ring(I)) && typeof(coefficient_ring(base_ring(I))) == Nemo.GaloisField
		standard_basis(I, ordering=start_ordering, complete_reduction=true, algorithm=:f4)
	else
		standard_basis(I, ordering=start_ordering, complete_reduction=true)
	end
	@assert dim(I) == 0
	I.gb[destination_ordering] = fglm(I.gb[start_ordering], ordering=destination_ordering)
end

