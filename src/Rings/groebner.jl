# groebner stuff #######################################################
@doc raw"""
    groebner_assure(I::MPolyIdeal, complete_reduction::Bool = false, need_global::Bool = false)
    groebner_assure(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an ideal `I` in a multivariate polynomial ring this function assures that a
Gröbner basis w.r.t. the given monomial ordering is attached to `I` in `I.gb`.
It *currently* also ensures that the basis is defined on the Singular side in
`I.gb.S`, but this should not be relied upon: use `singular_assure(I.gb)` before
accessing `I.gb.S`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

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
function groebner_assure(I::MPolyIdeal, complete_reduction::Bool = false, need_global::Bool = false)
	if !isempty(I.gb)
		for G in values(I.gb)
			need_global || return G
			is_global(G.ord) || continue
			complete_reduction || return G
			if !G.isReduced
				I.gb[G.ord] = _compute_standard_basis(G, G.ord, true)
				return I.gb[G.ord]
			end
		end
	end
	ord = default_ordering(base_ring(I))
	(need_global <= is_global(ord)) || error("Monomial ordering must be global.")
	I.gb[ord] = groebner_assure(I, ord, complete_reduction)
	return I.gb[ord]
end

function groebner_assure(I::MPolyIdeal, ordering::MonomialOrdering, complete_reduction::Bool = false)
    return get!(I.gb, ordering) do
        _compute_standard_basis(I.gens, ordering, complete_reduction)
    end
end

@doc raw"""
    _compute_standard_basis(B::IdealGens; ordering::MonomialOrdering,
                            complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `IdealGens` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a Gröbner basis (if `complete_reduction = true` the reduced Gröbner basis) of the
ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`. The Gröbner basis is then
returned in `B.S`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

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
@doc raw"""
    standard_basis(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool = false, algorithm::Symbol = :buchberger,
      weights::Vector{E} = ones(ngens(base_ring(I)))) where {E <: Integer}

Return a standard basis of `I` with respect to `ordering`.

The keyword `algorithm` can be set to
- `:buchberger` (implementation of Buchberger's algorithm in *Singular*),
- `:hilbert` (implementation of a Hilbert driven Gröbner basis computation in *Singular*),
- `:fglm` (implementation of the FGLM algorithm in *Singular*), and
- `:f4` (implementation of Faugère's F4 algorithm in the *msolve* package).

!!! note
    See the description of the functions `groebner_basis_hilbert_driven`, `fglm`, 
    and `f4` in the OSCAR documentation for some more details and for restrictions    
    on the input data when using these versions of the standard basis algorithm.

!!! note
    The returned standard basis is reduced if `ordering` is `global` and `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, ["x","y"]);

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
                        complete_reduction::Bool = false, algorithm::Symbol = :buchberger) 
  complete_reduction && @assert is_global(ordering)
  if haskey(I.gb, ordering) && (complete_reduction == false || I.gb[ordering].isReduced == true)
    return I.gb[ordering]
  end
  if algorithm == :buchberger
    if !haskey(I.gb, ordering)
      I.gb[ordering] = _compute_standard_basis(I.gens, ordering, complete_reduction)
    elseif complete_reduction == true
      I.gb[ordering] = _compute_standard_basis(I.gb[ordering], ordering, complete_reduction)
    end
  elseif algorithm == :fglm
    _compute_groebner_basis_using_fglm(I, ordering)
  elseif algorithm == :hilbert
    weights = _find_weights(gens(I))
    if !any(iszero, weights)
      J, target_ordering, hn = I, ordering, nothing
    else
      R = base_ring(I)
      K = iszero(characteristic(R)) && !haskey(I.gb, degrevlex(R)) ? _mod_rand_prime(I) : I
      S = base_ring(K)
      gb = groebner_assure(K, degrevlex(S))
      K_hom = homogenization(K, "w")
      gb_hom = IdealGens((p -> homogenization(p, base_ring(K_hom))).(gens(gb)))
      gb_hom.isGB = true
      K_hom.gb[degrevlex(S)] = gb_hom
      singular_assure(K_hom.gb[degrevlex(S)])
      hn = hilbert_series(quo(base_ring(K_hom), K_hom)[1])[1]
      J = homogenization(I, "w")
      weights = ones(Int, ngens(base_ring(J)))
      target_ordering = _extend_mon_order(ordering, base_ring(J))
    end
    GB = groebner_basis_hilbert_driven(J, destination_ordering=target_ordering,
                                       complete_reduction=complete_reduction,
                                       weights=weights,
                                       hilbert_numerator=hn)
    if base_ring(I) == base_ring(J)
      I.gb[ordering] = GB
    else
      GB_dehom_gens = [dehomogenization(p, base_ring(I), 1) for p in gens(GB)]
      I.gb[ordering] = IdealGens(GB_dehom_gens, ordering, isGB = true)
    end
  elseif algorithm == :f4
    groebner_basis_f4(I, complete_reduction=complete_reduction)
  end
  return I.gb[ordering]
end

@doc raw"""
    groebner_basis(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool = false, algorithm::Symbol = :buchberger)

If `ordering` is global, return a Gröbner basis of `I` with respect to `ordering`.

The keyword `algorithm` can be set to
- `:buchberger` (implementation of Buchberger's algorithm in *Singular*),
- `:hilbert` (implementation of a Hilbert driven Gröbner basis computation in *Singular*),
- `:fglm` (implementation of the FGLM algorithm in *Singular*), and
- `:f4` (implementation of Faugère's F4 algorithm in the *msolve* package).

!!! note
    See the description of the functions `groebner_basis_hilbert_driven`, `fglm`, 
    and `f4` in the OSCAR documentation for some more details and for restrictions    
    on the input data when using these versions of the standard basis algorithm.

!!! note
    The returned Gröbner basis is reduced if `complete_reduction = true`.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> G = groebner_basis(I)
Gröbner basis with elements
1 -> y^2 - x*z
2 -> x*y - z
3 -> x^2 - y
with respect to the ordering
degrevlex([x, y, z])

julia> elements(G)
3-element Vector{QQMPolyRingElem}:
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
julia> R, (x, y) = graded_polynomial_ring(QQ, ["x", "y"], [1, 3]);

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

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

julia> V = [3*x^3*y+x^3+x*y^3+y^2*z^2, 2*x^3*z-x*y-x*z^3-y^4-z^2,
               2*x^2*y*z-2*x*y^2+x*z^2-y^4];

julia> I = ideal(R, V);

julia> G = groebner_basis(I, ordering = lex(R), algorithm = :fglm);

julia> length(G)
8

julia> total_degree(G[8])
34

julia> leading_coefficient(G[8])
-91230304237130414552564280286681870842473427917231798336639893796481988733936505735341479640589040146625319419037353645834346047404145021391726185993823650399589880820226804328750
```
"""
function groebner_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)), complete_reduction::Bool=false,
                        algorithm::Symbol = :buchberger)
    is_global(ordering) || error("Ordering must be global")
    return standard_basis(I, ordering=ordering, complete_reduction=complete_reduction, algorithm=algorithm)
end

@doc raw"""
	groebner_basis_f4(I::MPolyIdeal, <keyword arguments>)

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
julia> R,(x,y,z) = polynomial_ring(GF(101), ["x","y","z"], ordering=:degrevlex)
(Multivariate polynomial ring in 3 variables over GF(101), fpMPolyRingElem[x, y, z])

julia> I = ideal(R, [x+2*y+2*z-1, x^2+2*y^2+2*z^2-x, 2*x*y+2*y*z-y])
ideal(x + 2*y + 2*z + 100, x^2 + 2*y^2 + 2*z^2 + 100*x, 2*x*y + 2*y*z + 100*y)

julia> groebner_basis_f4(I)
Gröbner basis with elements
1 -> x + 2*y + 2*z + 100
2 -> y*z + 82*z^2 + 10*y + 40*z
3 -> y^2 + 60*z^2 + 20*y + 81*z
4 -> z^3 + 28*z^2 + 64*y + 13*z
with respect to the ordering
degrevlex([x, y, z])
```
"""
function groebner_basis_f4(
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

@doc raw"""
    _compute_standard_basis_with_transform(B::BiPolyArray, ordering::MonomialOrdering, complete_reduction::Bool = false)

**Note**: Internal function, subject to change, do not use.

Given an `IdealGens` `B` and optional parameters `ordering` for a monomial ordering and `complete_reduction`
this function computes a standard basis (if `ordering` is a global monomial ordering and `complete_reduction = true`
the reduced Gröbner basis) of the ideal spanned by the elements in `B` w.r.t. the given monomial ordering `ordering`
and the transformation matrix from the ideal to the standard basis. Return value is a IdealGens together with a map.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

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

@doc raw"""
    standard_basis_with_transformation_matrix(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool=false)

Return a pair `G`, `T`, say, where `G` is a standard basis of `I` with respect to `ordering`, and `T` 
is a transformation matrix from `gens(I)` to `G`. That is, `gens(I)*T == G`.

!!! note
    The returned Gröbner basis is reduced if `ordering` is a global monomial odering and `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ,["x","y"]);

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

@doc raw"""
    groebner_basis_with_transformation_matrix(I::MPolyIdeal;
      ordering::MonomialOrdering = default_ordering(base_ring(I)),
      complete_reduction::Bool=false)

Return a pair `G`, `T`, say, where `G` is a Gröbner basis of `I` with respect to `ordering`, and `T` 
is a transformation matrix from `gens(I)` to `G`. That is, `gens(I)*T == G`.

!!! note
    The returned Gröbner basis is reduced if `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ,["x","y"]);

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
@doc raw"""
    syzygy_generators(G::Vector{<:MPolyRingElem})

Return generators for the syzygies on the polynomials given as elements of `G`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> S = syzygy_generators([x^3+y+2,x*y^2-13*x^2,y-14])
3-element Vector{FreeModElem{QQMPolyRingElem}}:
 (-y + 14)*e[2] + (-13*x^2 + x*y^2)*e[3]
 (-169*y + 2366)*e[1] + (-13*x*y + 182*x - 196*y + 2744)*e[2] + (13*x^2*y^2 - 2548*x^2 + 196*x*y^2 + 169*y + 338)*e[3]
 (-13*x^2 + 196*x)*e[1] + (-x^3 - 16)*e[2] + (x^4*y + 14*x^4 + 13*x^2 + 16*x*y + 28*x)*e[3]
```
"""
function syzygy_generators(a::Vector{<:MPolyRingElem})
  I = ideal(a)
  singular_assure(I)
  s = Singular.syz(I.gens.S)
  F = free_module(parent(a[1]), length(a))
  @assert rank(s) == length(a)
  return [F(s[i]) for i=1:Singular.ngens(s)]
end

# leading ideal #######################################################
@doc raw"""
    leading_ideal(G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(G[1]))) 
                                where T <: MPolyRingElem

Return the leading ideal of `G` with respect to `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5], ordering=degrevlex(R))
ideal(x*y^2, y^5)

julia> L = leading_ideal([x*y^2-3*x, x^3-14*y^5], ordering=lex(R))
ideal(x*y^2, x^3)
```
"""
function leading_ideal(G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(G[1]))) where { T <: MPolyRingElem }
    return ideal(parent(G[1]), [leading_monomial(f; ordering = ordering) for f in G])
end

function leading_ideal(I::IdealGens{T}) where { T <: MPolyRingElem }
    return ideal(base_ring(I), [leading_monomial(f; ordering = I.ord) for f in I])
end

function leading_ideal(I::IdealGens{T}, ordering::MonomialOrdering) where T <: MPolyRingElem
    return ideal(base_ring(I), [leading_monomial(f; ordering = ordering) for f in I])
end


@doc raw"""
    leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)))

Return the leading ideal of `I` with respect to `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

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

@doc raw"""
    normal_form_internal(I::Singular.sideal, J::MPolyIdeal, o::MonomialOrdering)

**Note**: Internal function, subject to change, do not use.

Compute the normal form of the generators `gens(I)` of the ideal `I` w.r.t. a
Gröbner basis of `J` and the monomial ordering `o`.

CAVEAT: This computation needs a Gröbner basis of `J` and the monomial ordering
`o`. If this Gröbner basis is not available, one is computed automatically.
This may take some time.

# Examples
```jldoctest
julia> R,(a,b,c) = polynomial_ring(QQ,["a","b","c"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> gens(groebner_basis(J))
2-element Vector{QQMPolyRingElem}:
 b + c - 1
 a*c - 2*a + c

julia> SR = singular_poly_ring(base_ring(J))
Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C)

julia> I = Singular.Ideal(SR,[SR(-1+c+b+a^3),SR(-1+b+c*a+2*a^3),SR(5+c*b+c^2*a)])
Singular ideal over Singular Polynomial Ring (QQ),(a,b,c),(dp(3),C) with generators (a^3 + b + c - 1, 2*a^3 + a*c + b - 1, a*c^2 + b*c + 5)

julia> Oscar.normal_form_internal(I,J,default_ordering(base_ring(J)))
3-element Vector{QQMPolyRingElem}:
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

@doc raw"""
	reduce(I::IdealGens, J::IdealGens; 
          ordering::MonomialOrdering = default_ordering(base_ring(J)))

Return a `Vector` whose elements are the underlying elements of `I`
reduced by the underlying generators of `J` w.r.t. the monomial
ordering `ordering`. `J` need not be a Gröbner basis. The returned
`Vector` will have the same number of elements as `I`, even if they
are zero.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(GF(11), ["x", "y", "z"]);

julia> I = ideal(R, [x^2, x*y - y^2]);

julia> J = ideal(R, [y^3])
ideal(y^3)

julia> reduce(J.gens, I.gens)
1-element Vector{fpMPolyRingElem}:
 y^3

julia> reduce(J.gens, groebner_basis(I))
1-element Vector{fpMPolyRingElem}:
 0

julia> reduce(y^3, [x^2, x*y-y^3])
x*y

julia> reduce(y^3, [x^2, x*y-y^3], ordering=lex(R))
y^3

julia> reduce([y^3], [x^2, x*y-y^3], ordering=lex(R))
1-element Vector{fpMPolyRingElem}:
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

@doc raw"""
	reduce(g::T, F::Vector{T}; 
           ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem

If `ordering` is global, return the remainder in a standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.
Otherwise, return the remainder in a *weak* standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.

	reduce(G::Vector{T}, F::Vector{T};
           ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem

Return a `Vector` which contains, for each element `g` of `G`, a remainder as above.

!!! note
    In the global case, the returned remainders are fully reduced.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"]);

julia> reduce(y^3, [x^2, x*y-y^3])
x*y

julia> reduce(y^3, [x^2, x*y-y^3], ordering = lex(R))
y^3
```

```jldoctest
julia> R, (z, y, x) = polynomial_ring(QQ, ["z", "y", "x"]);

julia> f1 = y-x^2; f2 = z-x^3;

julia> g = x^3*y-3*y^2*z^2+x*y*z;

julia> reduce(g, [f1, f2], ordering = lex(R))
-3*x^10 + x^6 + x^5
```
"""
function reduce(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(f))) where {T <: MPolyRingElem}
	@assert parent(f) == parent(F[1])
	R = parent(f)
	I = IdealGens(R, [f], ordering)
	J = IdealGens(R, F, ordering)
	redv = reduce(I, J, ordering=ordering)
	return redv[1]
end

function reduce(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyRingElem}
	@assert parent(F[1]) == parent(G[1])
	R = parent(F[1])
	I = IdealGens(R, F, ordering)
	J = IdealGens(R, G, ordering)
	return reduce(I, J, ordering=ordering)
end

@doc raw"""
	reduce_with_quotients_and_unit(g::T, F::Vector{T};
           ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem

Return the unit, the quotients and the remainder in a weak standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.

	reduce_with_quotients_and_unit(G::Vector{T}, F::Vector{T};
           ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem

Return a `Vector` which contains, for each element `g` of `G`, a unit, quotients, and a remainder as above.

!!! note
    In the global case, a standard representation with a fully reduced remainder is computed.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, ["x", "y", "z"]);

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
2-element Vector{QQMPolyRingElem}:
 -z^9 + z^7 + z^6 + z^4
 -3*z^7 + z^6

julia> U*G == Q*[f1, f2, f3]+H
true
```
"""
function reduce_with_quotients_and_unit(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyRingElem}
	@assert parent(f) == parent(F[1])
	R = parent(f)
	I = IdealGens(R, [f], ordering)
	J = IdealGens(R, F, ordering)
	u, q, r = _reduce_with_quotients_and_unit(I, J, ordering)
	return u, q, r[1]
end

function reduce_with_quotients_and_unit(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyRingElem}
	@assert parent(F[1]) == parent(G[1])
	R = parent(F[1])
	I = IdealGens(R, F, ordering)
	J = IdealGens(R, G, ordering)
	return _reduce_with_quotients_and_unit(I, J, ordering)
end

@doc raw"""
        reduce_with_quotients_and_unit(I::IdealGens, J::IdealGens; 
          ordering::MonomialOrdering = default_ordering(base_ring(J)))

Return a `Tuple` consisting of a `Generic.MatSpaceElem` `M`, a
`Vector` `res` whose elements are the underlying elements of `I`
reduced by the underlying generators of `J` w.r.t. the monomial
ordering `ordering` and a diagonal matrix `units` such that `M *
gens(J) + res == units * gens(I)`. If `ordering` is global then
`units` will always be the identity matrix, see also
`reduce_with_quotients`. `J` need not be a Gröbner basis. `res` will
have the same number of elements as `I`, even if they are zero.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(GF(11), ["x", "y"]);

julia> I = ideal(R, [x]);

julia> R, (x, y) = polynomial_ring(GF(11), ["x", "y"]);

julia> I = ideal(R, [x]);

julia> J = ideal(R, [x+1]);

julia> unit, M, res = reduce_with_quotients_and_unit(I.gens, J.gens, ordering = neglex(R))
([x+1], [x], fpMPolyRingElem[0])

julia> M * gens(J) + res == unit * gens(I)
true

julia> f = x^3*y^2-y^4-10
x^3*y^2 + 10*y^4 + 1

julia> F = [x^2*y-y^3, x^3-y^4]
2-element Vector{fpMPolyRingElem}:
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


@doc raw"""
        reduce_with_quotients(I::IdealGens, J::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(J)))

Return a `Tuple` consisting of a `Generic.MatSpaceElem` `M` and a
`Vector` `res` whose elements are the underlying elements of `I`
reduced by the underlying generators of `J` w.r.t. the monomial
ordering `ordering` such that `M * gens(J) + res == gens(I)` if `ordering` is global.
If `ordering` is local then this equality holds after `gens(I)` has been multiplied
with an unknown diagonal matrix of units, see reduce_with_quotients_and_unit` to
obtain this matrix. `J` need not be a Gröbner basis. `res` will have the same number
of elements as `I`, even if they are zero.

# Examples
```jldoctest
julia> R, (x, y, z) = polynomial_ring(GF(11), ["x", "y", "z"]);

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
([1 0 0; 0 0 1], fpMPolyRingElem[y^2, 0])

julia> M * gens(gb) + res == gens(I)
true

julia> f = x^3*y^2-y^4-10
x^3*y^2 + 10*y^4 + 1

julia> F = [x^2*y-y^3, x^3-y^4]
2-element Vector{fpMPolyRingElem}:
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

@doc raw"""
	reduce_with_quotients(g::T, F::Vector{T}; 
           ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem

If `ordering` is global, return the quotients and the remainder in a standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.
Otherwise, return the quotients and the remainder in a *weak* standard representation for `g` on division by the polynomials in `F` with respect to `ordering`.

	reduce_with_quotients(G::Vector{T}, F::Vector{T}; 
           ordering::MonomialOrdering = default_ordering(parent(F[1]))) where T <: MPolyRingElem

Return a `Vector` which contains, for each element `g` of `G`, quotients and a remainder as above.

!!! note
    In the global case, the returned remainders are fully reduced.

# Examples

```jldoctest
julia> R, (z, y, x) = polynomial_ring(QQ, ["z", "y", "x"]);

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
2-element Vector{QQMPolyRingElem}:
 -3*x^10 + x^6 + x^5
 -3*x^12 + x^7

julia> G == Q*[f1, f2]+H
true
```
"""
function reduce_with_quotients(f::T, F::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyRingElem}
	@assert parent(f) == parent(F[1])
	R = parent(f)
	I = IdealGens(R, [f], ordering)
	J = IdealGens(R, F, ordering)
	_, q, r = _reduce_with_quotients_and_unit(I, J, ordering)
	return q, r[1]
end

function reduce_with_quotients(F::Vector{T}, G::Vector{T}; ordering::MonomialOrdering = default_ordering(parent(F[1]))) where {T <: MPolyRingElem}
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

@doc raw"""
    normal_form(g::T, I::MPolyIdeal; 
      ordering::MonomialOrdering = default_ordering(base_ring(I))) where T <: MPolyRingElem

Compute the normal form of `g` mod `I` with respect to `ordering`.

    normal_form(G::Vector{T}, I::MPolyIdeal; 
      ordering::MonomialOrdering = default_ordering(base_ring(I))) where T <: MPolyRingElem

Return a `Vector` which contains for each element `g` of `G` a normal form as above.

# Examples
```jldoctest
julia> R,(a,b,c) = polynomial_ring(QQ,["a","b","c"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[a, b, c])

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> gens(groebner_basis(J))
2-element Vector{QQMPolyRingElem}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(-1+c+b+a^3, J)
a^3

julia> R,(a,b,c) = polynomial_ring(QQ,["a","b","c"])
(Multivariate polynomial ring in 3 variables over QQ, QQMPolyRingElem[a, b, c])

julia> A = [-1+c+b+a^3,-1+b+c*a+2*a^3,5+c*b+c^2*a]
3-element Vector{QQMPolyRingElem}:
 a^3 + b + c - 1
 2*a^3 + a*c + b - 1
 a*c^2 + b*c + 5

julia> J = ideal(R,[-1+c+b,-1+b+c*a+2*a*b])
ideal(b + c - 1, 2*a*b + a*c + b - 1)

julia> gens(groebner_basis(J))
2-element Vector{QQMPolyRingElem}:
 b + c - 1
 a*c - 2*a + c

julia> normal_form(A, J)
3-element Vector{QQMPolyRingElem}:
 a^3
 2*a^3 + 2*a - 2*c
 4*a - 2*c^2 - c + 5
```
"""
function normal_form(f::T, J::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(J))) where { T <: MPolyRingElem }
    singular_assure(J, ordering)
    I = Singular.Ideal(J.gens.Sx, J.gens.Sx(f))
    N = normal_form_internal(I, J, ordering)
    return N[1]
end

function normal_form(A::Vector{T}, J::MPolyIdeal; ordering::MonomialOrdering=default_ordering(base_ring(J))) where { T <: MPolyRingElem }
    singular_assure(J, ordering)
    I = Singular.Ideal(J.gens.Sx, [J.gens.Sx(x) for x in A])
    normal_form_internal(I, J, ordering)
end

@doc raw"""
    is_standard_basis(F::IdealGens; ordering::MonomialOrdering=default_ordering(base_ring(F)))

Tests if a given IdealGens `F` is a standard basis w.r.t. the given monomial ordering `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R,[x^2+y,x*y-y])
ideal(x^2 + y, x*y - y)

julia> is_standard_basis(I.gens, ordering=neglex(R))
false

julia> standard_basis(I, ordering=neglex(R))
Standard basis with elements
1 -> y
2 -> x^2
with respect to the ordering
neglex([x, y])

julia> is_standard_basis(I.gb[neglex(R)], ordering=neglex(R))
true
```
"""
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

@doc raw"""
    is_groebner_basis(F::IdealGens; ordering::MonomialOrdering=default_ordering(base_ring(F)))

Tests if a given IdealGens `F` is a Gröbner basis w.r.t. the given monomial ordering `ordering`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R,[x^2+y,x*y-y])
ideal(x^2 + y, x*y - y)

julia> is_groebner_basis(I.gens, ordering=lex(R))
false

julia> groebner_basis(I, ordering=lex(R))
Gröbner basis with elements
1 -> y^2 + y
2 -> x*y - y
3 -> x^2 + y
with respect to the ordering
lex([x, y])

julia> is_groebner_basis(I.gb[lex(R)], ordering=lex(R))
true
```
"""
function is_groebner_basis(F::IdealGens; ordering::MonomialOrdering = default_ordering(base_ring(F)))
    is_global(ordering) || error("Ordering must be global")
	return is_standard_basis(F, ordering=ordering)
end

@doc raw"""
    _fglm(G::IdealGens; ordering::MonomialOrdering)

Converts a Gröbner basis `G` w.r.t. a given global monomial ordering for `<G>`
to a Gröbner basis `H` w.r.t. another monomial ordering `ordering` for `<G>`.

**NOTE**: `_fglm` assumes that `G` is a reduced Gröbner basis (i.e. w.r.t. a global monomial ordering) and that `ordering` is a global monomial ordering.

# Examples
```jldoctest
julia> R, (x1, x2, x3, x4) = polynomial_ring(GF(101), ["x1", "x2", "x3", "x4"])
(Multivariate polynomial ring in 4 variables over GF(101), fpMPolyRingElem[x1, x2, x3, x4])

julia> J = ideal(R, [x1+2*x2+2*x3+2*x4-1,
       x1^2+2*x2^2+2*x3^2+2*x4^2-x1,
       2*x1*x2+2*x2*x3+2*x3*x4-x2,
       x2^2+2*x1*x3+2*x2*x4-x3
       ])
ideal(x1 + 2*x2 + 2*x3 + 2*x4 + 100, x1^2 + 100*x1 + 2*x2^2 + 2*x3^2 + 2*x4^2, 2*x1*x2 + 2*x2*x3 + 100*x2 + 2*x3*x4, 2*x1*x3 + x2^2 + 2*x2*x4 + 100*x3)

julia> groebner_basis(J, ordering=degrevlex(R), complete_reduction=true)
Gröbner basis with elements
1 -> x1 + 2*x2 + 2*x3 + 2*x4 + 100
2 -> x3^2 + 2*x2*x4 + 19*x3*x4 + 76*x4^2 + 72*x2 + 86*x3 + 42*x4
3 -> x2*x3 + 99*x2*x4 + 40*x3*x4 + 11*x4^2 + 65*x2 + 58*x3 + 30*x4
4 -> x2^2 + 2*x2*x4 + 30*x3*x4 + 45*x4^2 + 43*x2 + 72*x3 + 86*x4
5 -> x3*x4^2 + 46*x4^3 + 28*x2*x4 + 16*x3*x4 + 7*x4^2 + 58*x2 + 63*x3 + 15*x4
6 -> x2*x4^2 + 67*x4^3 + 56*x2*x4 + 58*x3*x4 + 45*x4^2 + 14*x2 + 86*x3
7 -> x4^4 + 65*x4^3 + 26*x2*x4 + 47*x3*x4 + 71*x4^2 + 37*x2 + 79*x3 + 100*x4
with respect to the ordering
degrevlex([x1, x2, x3, x4])

julia> Oscar._fglm(J.gb[degrevlex(R)], lex(R))
Gröbner basis with elements
1 -> x4^8 + 36*x4^7 + 95*x4^6 + 39*x4^5 + 74*x4^4 + 7*x4^3 + 45*x4^2 + 98*x4
2 -> x3 + 53*x4^7 + 93*x4^6 + 74*x4^5 + 26*x4^4 + 56*x4^3 + 15*x4^2 + 88*x4
3 -> x2 + 25*x4^7 + 57*x4^6 + 13*x4^5 + 16*x4^4 + 78*x4^3 + 31*x4^2 + 16*x4
4 -> x1 + 46*x4^7 + 3*x4^6 + 28*x4^5 + 17*x4^4 + 35*x4^3 + 9*x4^2 + 97*x4 + 100
with respect to the ordering
lex([x1, x2, x3, x4])
```
"""
function _fglm(G::IdealGens, ordering::MonomialOrdering)
	(G.isGB == true && G.isReduced == true) || error("Input must be a reduced Gröbner basis.") 
	singular_assure(G)
	Singular.dimension(G.S) == 0 || error("Dimension of corresponding ideal must be zero.")
	SR_destination, = Singular.polynomial_ring(base_ring(G.Sx),["$i" for i in gens(G.Sx)]; ordering = Singular.ordering_as_symbol(singular(ordering)))

	ptr = Singular.libSingular.fglmzero(G.S.ptr, G.Sx.ptr, SR_destination.ptr)
	return IdealGens(base_ring(G), Singular.sideal{Singular.spoly}(SR_destination, ptr, true))
end

@doc raw"""
    fglm(I::MPolyIdeal; start_ordering::MonomialOrdering = default_ordering(base_ring(I)),
                        destination_ordering::MonomialOrdering)

Given a **zero-dimensional** ideal `I`, return the reduced Gröbner basis of `I` with respect to `destination_ordering`.

!!! note
    Both `start_ordering` and `destination_ordering` must be global and the base ring of `I` must be a polynomial ring over a field.

!!! note
    The function implements the Gröbner basis conversion algorithm by **F**augère, **G**ianni, **L**azard, and **M**ora. See [FGLM93](@cite) for more information.

# Examples
```jldoctest
julia> R, (a, b, c, d, e) = polynomial_ring(QQ, ["a", "b", "c", "d", "e"]);

julia> f1 = a+b+c+d+e;

julia> f2 = a*b+b*c+c*d+a*e+d*e;

julia> f3 = a*b*c+b*c*d+a*b*e+a*d*e+c*d*e;

julia> f4 = b*c*d+a*b*c*e+a*b*d*e+a*c*d*e+b*c*d*e;

julia> f5 = a*b*c*d*e-1;

julia> I = ideal(R, [f1, f2, f3, f4, f5]);

julia> G = fglm(I, destination_ordering = lex(R));

julia> length(G)
8

julia> total_degree(G[8])
60

julia> leading_coefficient(G[8])
83369589588385815165248207597941242098312973356252482872580035860533111990678631297423089011608753348453253671406641805924218003925165995322989635503951507226650115539638517111445927746874479234
```
"""
function fglm(I::MPolyIdeal; start_ordering::MonomialOrdering = default_ordering(base_ring(I)), destination_ordering::MonomialOrdering)
	isa(coefficient_ring(I), AbstractAlgebra.Field) || error("The FGLM algorithm requires a coefficient ring that is a field.")
	(is_global(start_ordering) && is_global(destination_ordering)) || error("Start and destination orderings must be global.")
	haskey(I.gb, destination_ordering) && return I.gb[destination_ordering]
	if !haskey(I.gb, start_ordering)
		standard_basis(I, ordering=start_ordering, complete_reduction=true)
	elseif I.gb[start_ordering].isReduced == false
		I.gb[start_ordering] = _compute_standard_basis(I.gb[start_ordering], start_ordering, true)
	end

	I.gb[destination_ordering] = _fglm(I.gb[start_ordering], destination_ordering)

	return I.gb[destination_ordering]
end

@doc raw"""
    _compute_groebner_basis_using_fglm(I::MPolyIdeal, destination_ordering::MonomialOrdering)

Computes a reduced Gröbner basis for `I` w.r.t. `destination_ordering` using the FGLM algorithm.

**Note**: Internal function, subject to change, do not use.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, ["x", "y"])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> I = ideal(R,[x^2+y,x*y-y])
ideal(x^2 + y, x*y - y)

julia> Oscar._compute_groebner_basis_using_fglm(I, lex(R))
Gröbner basis with elements
1 -> y^2 + y
2 -> x*y - y
3 -> x^2 + y
with respect to the ordering
lex([x, y])

julia> I.gb[lex(R)]
Gröbner basis with elements
1 -> y^2 + y
2 -> x*y - y
3 -> x^2 + y
with respect to the ordering
lex([x, y])

julia> I.gb[degrevlex(R)]
Gröbner basis with elements
1 -> y^2 + y
2 -> x*y - y
3 -> x^2 + y
with respect to the ordering
degrevlex([x, y])
```
"""
function _compute_groebner_basis_using_fglm(I::MPolyIdeal,
	destination_ordering::MonomialOrdering)
	isa(coefficient_ring(I), AbstractAlgebra.Field) || error("The FGLM algorithm requires a coefficient ring that is a field.")
	haskey(I.gb, destination_ordering) && return I.gb[destination_ordering]
	is_global(destination_ordering) || error("Destination ordering must be global.")
	G = groebner_assure(I, true, true)
	start_ordering = G.ord
	dim(I) == 0 || error("Dimension of ideal must be zero.")
	I.gb[destination_ordering] = _fglm(G, destination_ordering)
end

@doc raw"""
    groebner_basis_hilbert_driven(I::MPolyIdeal{P}; destination_ordering::MonomialOrdering,
                    complete_reduction::Bool = false,
                    weights::Vector{Int} = ones(Int, ngens(base_ring(I))),
                    hilbert_numerator::Union{Nothing, ZZPolyRingElem} = nothing) 
                    where {P <: MPolyRingElem}

Return a Gröbner basis of `I` with respect to `destination_ordering`.

!!! note
    The function implements a version of the Hilbert driven Gröbner basis algorithm.
    See the corresponding section of the OSCAR documentation for some details.

!!! note
    All weights must be positive. If no weight vector is entered by the user, all weights 
    are set to 1. An error is thrown if the generators of `I` are not homogeneous with 
    respect to the corresponding (weighted) degree.  

!!! note
    If $R$ denotes the parent ring of $I$, and $p, q\in\mathbb Z[t]$ are polynomials
    such that $p/q$ represents the Hilbert series of $R/I$ as a rational function with 
    denominator $q = (1-t^{w_1})\cdots (1-t^{w_n}),$ where $n$ is the number of variables 
    of $R$, and $w_1, \dots, w_n$ are the assigned weights, then `hilbert_numerator` is 
    meant to be $p$. If this numerator is not entered by the user, it will be computed 
    internally.

# Examples
```jldoctest
julia> R, (a, b, c, d, e, f, g) = polynomial_ring(QQ, ["a", "b", "c", "d", "e", "f", "g"]);

julia> V = [-3*a^2+2*f*b+3*f*d, (3*g*b+3*g*e)*a-3*f*c*b,
                      -3*g^2*a^2-c*b^2*a-g^2*f*e-g^4, e*a-f*b-d*c];

julia> I = ideal(R, V);

julia> o = degrevlex([a, b, c])*degrevlex([d, e, f, g]);

julia> G = groebner_basis_hilbert_driven(I, destination_ordering = o);

julia> length(G)
296

julia> total_degree(G[49])
30
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(GF(32003), ["x", "y", "z"]);

julia> f1 = x^2*y+169*y^21+151*x*y*z^10;

julia> f2 = 6*x^2*y^4+x*z^14+3*z^24;

julia> f3 = 11*x^3+5*x*y^10*z^10+2*y^20*z^10+y^10*z^20;

julia> I = ideal(R, [f1, f2,f3]);

julia> W = [10, 1, 1];

julia> GB = groebner_basis_hilbert_driven(I, destination_ordering = lex(R), weights = W);

julia> length(GB)
40
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(GF(32003), ["x", "y", "z"]);

julia> f1 = x^2*y+169*y^21+151*x*y*z^10;

julia> f2 = 6*x^2*y^4+x*z^14+3*z^24;

julia> f3 = 11*x^3+5*x*y^10*z^10+2*y^20*z^10+y^10*z^20;

julia> I = ideal(R, [f1, f2,f3]);

julia> W = [10, 1, 1];

julia> S, t = polynomial_ring(ZZ, "t")
(Univariate polynomial ring in t over ZZ, t)

julia> hn = -t^75 + t^54 + t^51 + t^45 - t^30 - t^24 - t^21 + 1
-t^75 + t^54 + t^51 + t^45 - t^30 - t^24 - t^21 + 1

julia> GB = groebner_basis_hilbert_driven(I, destination_ordering = lex(R), weights = W, hilbert_numerator = hn);

julia> length(GB)
40
```
"""
function groebner_basis_hilbert_driven(I::MPolyIdeal{P};
                                       destination_ordering::MonomialOrdering,
                                       complete_reduction::Bool = false,
                                       weights::Vector{Int} = ones(Int, ngens(base_ring(I))),
                                       hilbert_numerator::Union{Nothing, ZZPolyRingElem} = nothing) where {P <: MPolyRingElem}
  
  all(f -> _is_homogeneous(f, weights), gens(I)) || error("I must be given by generators homogeneous with respect to the given weights.")
  isa(coefficient_ring(I), AbstractAlgebra.Field) || error("The underlying coefficient ring of I must be a field.")
  ordering = destination_ordering
  is_global(ordering) || error("Destination ordering must be global.")
  haskey(I.gb, ordering) && return I.gb[ordering]
  if isnothing(hilbert_numerator)
    if isempty(I.gb)
      J = iszero(characteristic(base_ring(I))) ? _mod_rand_prime(I) : I
      G = groebner_assure(J, wdegrevlex(base_ring(J), weights))
    else
      G = groebner_assure(I)
    end

    if characteristic(base_ring(I)) > 0 && ordering == wdegrevlex(base_ring(I), weights)
      return G
    end
    singular_assure(G)
    h = Singular.hilbert_series(G.S, weights)
  else
    # Quoting from the documentation of Singular.hilbert_series:
    # The coefficient vector is returned as a `Vector{Int32}`, and the last element is not actually part of the coefficients of Q(t).
    # what?
    h = (Int32).([coeff(hilbert_numerator, i) for i in 0:degree(hilbert_numerator)+1])
  end

  singular_assure(I.gens, ordering)
  singular_ring = I.gens.Sx
  J  = Singular.Ideal(singular_ring, gens(I.gens.S)...)
  i  = Singular.std_hilbert(J, h, (Int32).(weights),
                            complete_reduction = complete_reduction)
  GB = IdealGens(I.gens.Ox, i, complete_reduction)
  GB.isGB = true
  GB.ord = ordering
  if isdefined(GB, :S)
    GB.S.isGB  = true
  end
  I.gb[destination_ordering] = GB
  return GB
end

# Helper functions for groebner_basis_with_hilbert

function _extract_weights(T::MPolyDecRing)
  if !is_z_graded(T)
    error("Ring must be graded by the Integers.")
  end
  return [Int(first(gr_elem.coeff)) for gr_elem in T.d]
end

function _extend_mon_order(ordering::MonomialOrdering,
                           homogenized_ring::MPolyDecRing)

  nvars = ngens(ordering.R)
  m = canonical_matrix(ordering)
  m_hom = similar(m, nvars + 1, nvars + 1)
  m_hom[1, :] = ones(Int, nvars + 1)
  m_hom[2:end, 2:end] = m
  return matrix_ordering(homogenized_ring, m_hom)
end

function _mod_rand_prime(I::MPolyIdeal)
  p = 32771
  while true
    p = Hecke.next_prime(p)
    
    base_field = GF(p)
    ModP, _ = polynomial_ring(base_field, ngens(base_ring(I)))
    I_mod_p_gens =
      try
        [map_coefficients(base_field, f; parent=ModP) for f in gens(I)]
      catch e
        # this precise error is thrown if the chosen prime p divides one
        # of the denominators of the coefficients of the generators of I.
        # In this case we simply choose the next prime and try again.
        if e == ErrorException("Unable to coerce") 
          continue
        else
          rethrow(e)
        end
      end
    return ideal(ModP, I_mod_p_gens)
  end
end

# check homogeneity w.r.t. some weights

function _is_homogeneous(f::MPolyRingElem, weights::Vector{Int})
  w = sum(weights .* first(exponents(f)))
  all(sum(weights .* e) == w for e in exponents(f))
end


# check homogeneity w.r.t. total degree
function _is_homogeneous(f::MPolyRingElem)
  leadexpv,tailexpvs = Iterators.peel(AbstractAlgebra.exponent_vectors(f))
  d = sum(leadexpv)
  for tailexpv in tailexpvs
    if d!=sum(tailexpv)
      return false
    end
  end
  return true
end
  

# compute weights such that F is a homogeneous system w.r.t. these weights
function _find_weights(F::Vector{P}) where {P <: MPolyRingElem}

  if all(_is_homogeneous, F)
    return ones(Int, ngens(parent(F[1])))
  end

  nrows = sum((length).(F)) - length(F)
  ncols = ngens(parent(first(F)))

  exp_diffs = permutedims(reduce(hcat, [e[i] - e[1] for e in
                                          (collect).((exponents).(F))
                                          for i in 2:length(e)]))
  K = kernel(matrix(QQ, nrows, ncols, exp_diffs))[2]
  isempty(K) && return zeros(Int, ncols)
  # Here we try to find a vector with strictly positive entries in K
  # this method to find such a vector is taken from
  # https://mathoverflow.net/questions/363181/intersection-of-a-vector-subspace-with-a-cone
  Pol = polyhedron(-K,  zeros(Int, ncols))
  !is_feasible(Pol) && return zeros(Int, ncols)
  pos_vec = zeros(Int, ncols)
  for i in 1:ncols
    ei = [j == i ? one(QQ) : zero(QQ) for j in 1:ncols]
    obj_func = ei * K
    L = linear_program(Pol, obj_func)
    m, v = solve_lp(L)
    if isnothing(v)
      Pol_new = intersect(Pol, polyhedron(ei*K, [1]))
      L = linear_program(Pol_new, obj_func)
      v = optimal_vertex(L)
    end
    pos_vec += K*(v.p)
  end
  ret = (Int).(lcm((denominator).(pos_vec)) .* pos_vec)
  ret = (x -> div(x, gcd(ret))).(ret) 
  # assure that the weights fit in Int32 for singular
  return all(ret .< 2^32) ? ret : zeros(Int,ncols)
end
