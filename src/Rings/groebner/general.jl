# groebner stuff #######################################################
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
julia> R,(x,y) = polynomial_ring(QQ, [:x,:y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> A = Oscar.IdealGens([x*y-3*x,y^3-2*x^2*y])
Ideal generating system with elements
  1: x*y - 3*x
  2: -2*x^2*y + y^3

julia> B = Oscar._compute_standard_basis(A, degrevlex(R))
Gröbner basis with elements
  1: x*y - 3*x
  2: y^3 - 6*x^2
  3: 2*x^3 - 9*x
with respect to the ordering
  degrevlex([x, y])
```
"""
function _compute_standard_basis(B::IdealGens, ordering::MonomialOrdering, complete_reduction::Bool = false)
  gensSord = singular_generators(B, ordering)
  i = Singular.std(gensSord, complete_reduction = complete_reduction)
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
    standard_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)),
                   complete_reduction::Bool = false, algorithm::Symbol = :buchberger)

Return a standard basis of `I` with respect to `ordering`.

The keyword `algorithm` can be set to
- `:buchberger` (implementation of Buchberger's algorithm in *Singular*),
- `:modular` (implementation of multi-modular approach, if applicable),
- `:f4` (implementation of Faugère's F4 algorithm in the *msolve* package),
- `:fglm` (implementation of the FGLM algorithm in *Singular*),
- `:hc` (implementation of Buchberger's algorithm in *Singular* trying to first compute the highest corner modulo some prime), and
- `:hilbert` (implementation of a Hilbert driven Gröbner basis computation in *Singular*).
- `:markov` (implementation to compute Gröbner basis for lattice ideals in *4ti2*)

!!! note
    See the description of the functions `groebner_basis_hilbert_driven`, `fglm`, 
    and `f4` in the OSCAR documentation for some more details and for restrictions    
    on the input data when using these versions of the standard basis algorithm.

!!! note
    The returned standard basis is reduced if `ordering` is `global` and `complete_reduction = true`.

# Examples
```jldoctest
julia> R,(x,y) = polynomial_ring(QQ, [:x,:y]);

julia> I = ideal([x*(x+1), x^2-y^2+(x-2)*y]);

julia> standard_basis(I, ordering = negdegrevlex(R))
Standard basis with elements
  1: x
  2: y
with respect to the ordering
  negdegrevlex([x, y])
```
"""
function standard_basis(I::MPolyIdeal; ordering::MonomialOrdering = default_ordering(base_ring(I)),
                        complete_reduction::Bool = false, algorithm::Symbol = :buchberger)
  complete_reduction && @assert is_global(ordering)
  @req is_exact_type(elem_type(base_ring(I))) "This functionality is only supported over exact fields."
  if haskey(I.gb, ordering) && (complete_reduction == false || I.gb[ordering].isReduced == true)
    return I.gb[ordering]
  end
  if algorithm == :buchberger
    if !haskey(I.gb, ordering)
      I.gb[ordering] = _compute_standard_basis(I.gens, ordering, complete_reduction)
    elseif complete_reduction == true
      I.gb[ordering] = _compute_standard_basis(I.gb[ordering], ordering, complete_reduction)
    end
  elseif algorithm == :modular
    if is_f4_applicable(I, ordering)
      #  since msolve v0.7.0 is most of the time more efficient
      #  to compute a reduced GB by default
      groebner_basis_f4(I, complete_reduction=true)
    elseif base_ring(I) isa QQMPolyRing
      groebner_basis_modular_qq(I, ordering=ordering)
    elseif coefficient_ring(I) isa FracField{QQMPolyRingElem}
      groebner_basis_modular_qt(I, ordering=ordering)
      println(I.gb[ordering])
    else
      error("Modular option not applicable in this setting.")
    end
  elseif algorithm == :fglm
    _compute_groebner_basis_using_fglm(I, ordering)
  elseif algorithm == :hc
    standard_basis_highest_corner(I, ordering=ordering)
  elseif algorithm == :hilbert
    weights = _find_weights(gens(I))
    if !any(iszero, weights)
      J, target_ordering, hn = I, ordering, nothing
    else
      R = base_ring(I)
      K = iszero(characteristic(R)) && !haskey(I.gb, degrevlex(R)) ? _mod_rand_prime(I) : I
      S = base_ring(K)
      gb = groebner_assure(K, degrevlex(S))
      # 2024-02-09 Next lines "blindly" updated to use new homogenization UI
      H = homogenizer(S, "w")
      K_hom = H(K)
      gb_hom = IdealGens(H.(gens(gb)))
      gb_hom.isGB = true
      K_hom.gb[degrevlex(S)] = gb_hom
      singular_assure(K_hom.gb[degrevlex(S)])
      hn = hilbert_series(quo(base_ring(K_hom), K_hom)[1])[1]
      H2 = homogenizer(R, "w")
      J = H2(I)
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
      DH2 = dehomogenizer(H2)
      GB_dehom_gens = DH2.(gens(GB))
      I.gb[ordering] = IdealGens(GB_dehom_gens, ordering, isGB = true)
    end
  elseif algorithm == :f4
    #  since msolve v0.7.0 is most of the time more efficient
    #  to compute a reduced GB by default
    groebner_basis_f4(I, complete_reduction=true)
  elseif algorithm == :markov
    @req all(i -> length(i) == 2, gens(I)) "Ideal needs to be a lattice ideal to use markov algorithm"
    @req all(iszero, sum.(coefficients.(gens(I)))) "Ideal needs to be a lattice ideal to use markov algorithm"
    @req all(x -> all(is_unit.(x)), coefficients.(gens(I))) "Ideal needs to be a lattice ideal to use markov algorithm"
    I.gb[ordering] = _groebner4ti2(I, ordering)
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
- `:modular` (implementation of multi-modular approach, if applicable),
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
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [y-x^2, z-x^3]);

julia> G = groebner_basis(I)
Gröbner basis with elements
  1: y^2 - x*z
  2: x*y - z
  3: x^2 - y
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
  1: y^3 - z^2
  2: x*z - y^2
  3: x*y - z
  4: x^2 - y
with respect to the ordering
  lex([x, y, z])
```
```jldoctest
julia> R, (x, y) = graded_polynomial_ring(QQ, [:x, :y], [1, 3]);

julia> I = ideal(R, [x*y-3*x^4,y^3-2*x^6*y]);

julia> groebner_basis(I)
Gröbner basis with elements
  1: 3*x^4 - x*y
  2: 2*x^3*y^2 - 3*y^3
  3: x*y^3
  4: y^4
with respect to the ordering
  wdegrevlex([x, y], [1, 3])
```

```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

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

function is_known(
    ::typeof(groebner_basis), I::MPolyIdeal; 
    ordering::Union{MonomialOrdering, Nothing}=nothing
  )
  ordering === nothing && return !is_empty(I.gb) # Check whether there is any groebner basis
  return haskey(I.gb, ordering)
end

