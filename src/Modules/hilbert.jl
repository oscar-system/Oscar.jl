# This cannot be included in src/Rings/hilbert.jl because the include
# order defined in src/Oscar.jl includes ring stuff before module stuff
# (and it probably would not work if we inverted the inclusion order).

# ==================================================================
# Hilbert series numerator for modules

function _power_product(T, expv)
  return prod([T[k]^expv[k]  for k in 1:length(expv)]);
end


@doc raw"""
    hilbert_series(M::SubquoModule; parent::Union{Nothing,Ring} = nothing)

Compute a pair `(N,D)` where `N` and `D` are the non-reduced numerator and denominator of the Hilbert
series of the subquotient `M`.  If the kwarg `parent` is supplied `N` and `D` are computed in the ring `parent`.

!!! note
    Applied to a homogeneous subquotient `M`, the function first computes a Groebner basis to
    obtain the leading term module; the rest of the computation uses this latter module
    (sliced into ideals, one for each ambient free module component).

# Examples
```jldoctest
julia> R, _ = polynomial_ring(QQ, ["x", "y", "z"]);

julia> Z = abelian_group(0);

julia> Rg, (x, y, z) = grade(R, [Z[1],Z[1],Z[1]]);

julia> F = graded_free_module(Rg, 1);

julia> A = Rg[x; y];

julia> B = Rg[x^2; y^3; z^4];

julia> M = SubquoModule(F, A, B);

julia> hilbert_series(M)
((-t^9 + t^7 + 2*t^6 - t^5 - t^3 - 2*t^2 + 2*t, Factored element with data
Dict{AbstractAlgebra.Generic.LaurentMPolyWrap{ZZRingElem, ZZMPolyRingElem, AbstractAlgebra.Generic.LaurentMPolyWrapRing{ZZRingElem, ZZMPolyRing}}, ZZRingElem}(-t + 1 => 3)), (GrpAb: Z, Identity map with

Domain:
=======
GrpAb: Z))

```
"""
function HSNum_module(SubM::SubquoModule{T}, HSRing::Ring)  where T <: MPolyRingElem
  C,phi = present_as_cokernel(SubM, :with_morphism);  # phi is the morphism
  LM = leading_module(C.quo);
  F = ambient_free_module(C.quo);
  r = rank(F);
  P = base_ring(C.quo);
  # short-cut for module R^0 (otherwise defn if HSeriesRing below gives index error)
  if iszero(r)
    return HSNum_abbott(quo(P,ideal(P,[1]))[1], HSRing)
  end
  GensLM = gens(LM);
  L = [[] for _ in 1:r];  # L[k] is list of monomial gens for k-th cooord
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
###  # If parent === nothing, should we compute a HSNum for some ideal first then use its ring for all later calls?
  HSeriesList = [HSNum_abbott(quo(P,I)[1], HSRing)  for I in IdealList];
  shifts = [degree(phi(g))  for g in gens(F)];
  @vprintln :hilbert 1 "HSNum_module: shifts are $(shifts)";
  shift_expv = [gen_repr(d)  for d in shifts];
  @vprintln :hilbert 1 "HSNum_module: shift_expv are $(shift_expv)";
  HSeriesRing = Oscar.parent(HSeriesList[1]);
  @vprintln :hilbert 1 "HSNum_module: HSeriesRing = $(HSeriesRing)";
  t = gens(HSeriesRing);
  ScaleFactor = [_power_product(t,e)  for e in shift_expv];
  result = sum([ScaleFactor[k]*HSeriesList[k]  for k in 1:r]);
  return result;
end

# function HSNum_module(F::FreeMod{T}; parent::Union{Nothing,Ring} = nothing)  where T <: MPolyRingElem
#   # ASSUME F is graded free module
#   return HSNum_module(sub(F,gens(F))[1]; parent=parent)
# end

function hilbert_series(SubM::SubquoModule{T}; parent::Union{Nothing,Ring} = nothing #=, backend::Symbol = :Abbott=#)  where T <: MPolyRingElem
  R = base_ring(SubM)
  @req coefficient_ring(R) isa AbstractAlgebra.Field "The coefficient ring must be a field"
##???  @req is_positively_graded(R) "The base ring must be positively graded"

  # Wrap the case where G is abstractly isomorphic to ℤᵐ, but not realized as a 
  # free Abelian group. 
  #
  # We use the Smith normal form to get there, recreate the graded ring with the 
  # free grading group, do the computation there and return the isomorphism for 
  # the grading. 
  G = grading_group(R)
  if !is_zm_graded(R)
    H, iso = snf(G)
    V = [preimage(iso, x) for x in gens(G)]
    isoinv = hom(G, H, V)
#    W = R.d
    W = [isoinv(R.d[i]) for i = 1:length(W)]
    S, _ = graded_polynomial_ring(coefficient_ring(R), symbols(R), W)
    map_into_S = hom(R, S, gens(S))
    J = map_into_S(modulus(SubM))  ##?????
    SubM2, _ = quo(S, J)
#??    change_res = hom(A, AA, gens(AA); check=false)
#??    change_res_inv = hom(AA, A, gens(A); check=false)
    (numer, denom), _ = multi_hilbert_series(SubM2; algorithm, parent)
    return (numer, denom), (H, iso)
  end

  # Now we may assume that the grading group is free Abelian.
  m = ngens(G)  
  n = ngens(R)
  HSRing = _hilbert_series_ring(parent, m)
  # Get the weights as Int values: W[k] contain the weight(s) of x[k]
  W = [[ Int(R.d[i][j])  for j in 1:m]  for i in 1:n]
  fac_denom = _hilbert_series_denominator(HSRing, W)
  numer = HSNum_module(SubM, HSRing #=; backend=backend=#)
  return (numer, fac_denom), (G, identity_map(G))
end

function hilbert_series(F::FreeMod{T}; parent::Union{Nothing,Ring} = nothing)  where T <: MPolyRingElem
  return hilbert_series(sub(F,gens(F))[1]; parent=parent)
end
