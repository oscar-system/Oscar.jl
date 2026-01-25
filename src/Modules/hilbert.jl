# This cannot be included in src/Rings/hilbert.jl because the include
# order defined in src/Oscar.jl includes ring stuff before module stuff
# (and it probably would not work if we inverted the inclusion order).

# ==================================================================
# Hilbert series numerator for modules

function _power_product(T, expv)
  return prod([T[k]^expv[k]  for k in 1:length(expv)]);
end


function HSNum_module(SubM::SubquoModule{T}, HSRing::Ring, backend::Symbol=:Abbott)  where T <: MPolyRingElem
  C,phi = present_as_cokernel(SubM, :with_morphism);  # phi is the morphism
  LM = leading_module(C.quo);
  F = ambient_free_module(C.quo);
  r = rank(F);
  P = base_ring(C.quo);
  # short-cut for module R^0 (to avoid problems with empty sum below)
  if iszero(r)
    return multi_hilbert_series(quo(P,ideal(P, one(P)))[1]; parent=HSRing, backend=backend)[1][1]
  end
  GensLM = gens(LM);
  L = [[] for _ in 1:r];  # L[k] will be list of monomial gens for k-th coord
  # Nested loop below extracts the coordinate monomial ideals -- is there a better way?
  for g in GensLM
    SR = coordinates(ambient_representative(g)); # should have length = 1
    for j in 1:r
      if SR[j] != 0
        push!(L[j], SR[j]);
      end
    end
  end
  IdealList = [ideal(P,G)  for G in L];
  HSeriesList = [multi_hilbert_series(quo(P,I)[1]; parent=HSRing, backend=backend)[1][1]  for I in IdealList];
  shifts = [degree(phi(g))  for g in gens(F)];
  @vprintln :hilbert 1 "HSNum_module: shifts are $(shifts)";
  shift_expv = [gen_repr(d)  for d in shifts];
  @vprintln :hilbert 1 "HSNum_module: shift_expv are $(shift_expv)";
  t = gens(HSRing);
  ScaleFactor = [_power_product(t,e)  for e in shift_expv];
  result = sum([ScaleFactor[k]*HSeriesList[k]  for k in 1:r]);
  return result;
end


@doc raw"""
    multi_hilbert_series(M::SubquoModule; parent::Ring)

Compute a pair of pairs `(N ,D), (H ,iso)` where `N` and `D` are the non-reduced numerator and denominator of the Hilbert
series of the subquotient `M`, and `H` is the SNF of the grading group together with the identifying isomorphism `iso`.
If the kwarg `parent` is supplied `N` and `D` are computed in the ring `parent`.

!!! note
    CURRENT LIMITATION: the grading must be over ZZ^m (more general gradings are not yet supported)

!!! note
    Applied to a homogeneous subquotient `M`, the function first computes a Groebner basis to
    obtain the leading term module; the rest of the computation uses this latter module
    (sliced into ideals, one for each ambient free module component).

# Examples
```jldoctest
julia> Rg, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]; weights = [4,3,2]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> (num,den),_ = multi_hilbert_series(M);

julia> num
-t^25 + 2*t^17 + t^16 + t^15 - t^12 - t^11 - t^9 - t^8 - t^7 + t^4 + t^3

julia> den
(-t^4 + 1)^1*(-t^3 + 1)^1*(-t^2 + 1)^1

```
"""
function multi_hilbert_series(
    SubM::SubquoModule{T}; 
    parent::Ring = multi_hilbert_series_parent(base_ring(SubM)),
    backend::Symbol = :Abbott
  )  where T <: MPolyRingElem
  R = base_ring(SubM)
  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
  @req is_positively_graded(R) "ring must be positively graded"

  # Wrap the case where G is abstractly isomorphic to ℤᵐ, but not realized as a 
  # free Abelian group. 
  #
  # We use the Smith normal form to get there, recreate the graded ring with the 
  # free grading group, do the computation there and return the isomorphism for 
  # the grading. 
  G = grading_group(R)
  if !is_zm_graded(R)
    error("non-ZZ^m-grading not yet implemented for modules (see issue #2657)")  ### Needs improvement to change_base_ring (see below)
    H, iso = snf(G)
    V = [preimage(iso, x) for x in gens(G)]
    isoinv = hom(G, H, V)
    W = [isoinv(R.d[i]) for i = 1:ngens(R)]
    S, _ = graded_polynomial_ring(coefficient_ring(R), symbols(R); weights=W, cached=false)
    map_into_S = hom(R, S, gens(S))
    SubM2,_ = change_base_ring(map_into_S,SubM) # !!! BUG this seems to forget that things are graded BUG (issue #2657) !!!
    (numer, denom), _ = hilbert_series(SubM2; parent=parent, backend=backend)
    return (numer, denom), (H, iso)
  end

  # Now we may assume that the grading group is free Abelian.
  m = ngens(G)
  n = ngens(R)
  @req ngens(parent) >= m "Parent ring does not contain sufficiently many variables"
  # Get the weights as Int values: W[k] contains the weight(s) of x[k]
  W = [[ Int(R.d[i][j])  for j in 1:m]  for i in 1:n]
  denom = _hilbert_series_denominator(parent, W)
  numer = HSNum_module(SubM, parent, backend)
  return (numer, denom), (G, id_hom(G))
end

function multi_hilbert_series(
    F::FreeMod{T};
    parent::Ring = multi_hilbert_series_parent(base_ring(F)),
    backend::Symbol = :Abbott
  )  where T <: MPolyRingElem
  @req is_positively_graded(base_ring(F)) "ring must be positively graded"
  return multi_hilbert_series(sub_object(F,gens(F)); parent=parent, backend=backend)
end


@doc raw"""
    hilbert_series(M::SubquoModule; parent::Ring)

Compute a pair `(N,D)` where `N` and `D` are the non-reduced numerator and denominator of the Hilbert
series of the subquotient `M`.  If the kwarg `parent` is supplied `N` and `D` are computed in the ring `parent`.

!!! note
    Applied to a homogeneous subquotient `M`, the function first computes a Groebner basis to
    obtain the leading term module; the rest of the computation uses this latter module
    (sliced into ideals, one for each ambient free module component).

# Examples
```jldoctest
julia> R, _ = polynomial_ring(QQ, [:x, :y, :z]);

julia> Z = abelian_group(0);

julia> Rg, (x, y, z) = grade(R, [Z[1],Z[1],Z[1]]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> num,den = hilbert_series(M);

julia> num
-t^9 + t^7 + 2*t^6 - t^5 - t^3 - 2*t^2 + 2*t

julia> den
(-t + 1)^3

```
"""
function hilbert_series(
    SubM::SubquoModule{T};
    parent::Ring = hilbert_series_parent(base_ring(SubM)),
    backend::Symbol = :Abbott
  )  where T <: MPolyRingElem
  @req is_z_graded(base_ring(SubM)) "ring must be ZZ-graded; use `multi_hilbert_series` otherwise"
  HS, _ = multi_hilbert_series(SubM; parent=parent, backend=backend)
  return HS
end

function hilbert_series(
    F::FreeMod{T}; 
    parent::Ring = hilbert_series_parent(base_ring(F)), 
    backend::Symbol = :Abbott
  )  where T <: MPolyRingElem
  @req is_z_graded(base_ring(F)) "ring must be ZZ-graded; use `multi_hilbert_series` otherwise"
  return hilbert_series(sub_object(F,gens(F)); parent=parent, backend=backend)
end

function hilbert_series_parent(S::MPolyDecRing)
  if !isdefined(S, :hilbert_series_parent)
    S.hilbert_series_parent = laurent_polynomial_ring(ZZ, :t; cached=false)[1]
  end
  return S.hilbert_series_parent
end

function multi_hilbert_series_parent(S::MPolyDecRing)
  if !isdefined(S, :multi_hilbert_series_parent)
    G = grading_group(S)
    m = ngens(G)
    S.multi_hilbert_series_parent = laurent_polynomial_ring(ZZ, (isone(m) ? [:t] : (:t => 1:m)); cached=false)[1]
  end
  return S.multi_hilbert_series_parent
end
