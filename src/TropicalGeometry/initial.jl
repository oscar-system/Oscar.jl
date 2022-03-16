###
# Computing initial forms and initial ideals in Oscar
# ===================================================
###

@doc Markdown.doc"""
    valued_weighted_degree(f::MPolyElem, val::TropicalSemiringMap, w::Vector; pertubation::Vector=[], return_vector::Bool=false)

Return the valued weighted degree of a polynomial `f` with respect to valuation
`val` and weight vector `w`. In other words, returns the tropicalized
polynomial of `f` with respect to valuation `val` evaluated at `w`.

If `return_vector=true`, returns a vector whose i-th entry is the valued
weighted degree of the i-th term of `f`.

# Examples
```jldoctest
julia> Kxy, (x,y) = PolynomialRing(QQ,["x", "y"]);

julia> val_2 = TropicalSemiringMap(QQ,2);

julia> val_trivial = TropicalSemiringMap(QQ);

julia> w = [1,1];

julia> f = 2*x+2*y+1;

julia> valued_weighted_degree(f, val_2, w)
(0)

julia> valued_weighted_degree(f, val_2, w, return_vector=true)
((0), Oscar.TropicalSemiringElem{typeof(min)}[(2), (2), (0)])

julia> valued_weighted_degree(f, val_trivial, w, return_vector=true)
((0), Oscar.TropicalSemiringElem{typeof(min)}[(1), (1), (0)])

```
"""
function valued_weighted_degree(f::MPolyElem, val::TropicalSemiringMap, w::Vector; pertubation::Vector=[], return_vector::Bool=false)
  # compute the weighted degrees shifted by the coefficient valuations
  vwds = [val(c)*tropical_semiring(val)(dot(w,alpha)) for (c,alpha) in zip(coefficients(f),exponent_vectors(f))]

  # compute the minimal degree
  # (note: max tropical semiring is reversely ordered, so min in the semiring is max in the conventional sense)
  vwd = min(vwds...)

  if isempty(pertubation)
    # if no pertubation is specified, compute the maximum and return vector if required
    if return_vector
      return vwd,vwds
    end
    return vwd
  else
    # if pertubation is specified, then compute the pertubed degrees
    vwdsPerp = [tropical_semiring(val)(dot(pertubation,alpha)) for alpha in exponent_vectors(f)]
    # compute the mininum amongst all pertubations with maximal original degree
    vwdPerp = min([vwdsPerp[i] for i in 1:length(vwds) if vwds[i]==vwd]...)

    if return_vector
      return vwd,vwdPerp,vwds,vwdsPerp
    end
    return vwd,vwdPerp
  end

end
export valued_weighted_degree



# # not wrong, but not sure whether needed
# function weighted_degree(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, w::Vector; return_vector::Bool=false)
#   trivial_val = TropicalSemiringMap(coefficient_ring(f))
#   return valued_weighted_degree(f, trivial_val, w, return_vector=return_vector)
# end
# export weighted_degree



@doc Markdown.doc"""
    initial(f::MPolyElem, val::TropicalSemiringMap, w::Vector, convention::Union{typeof(min),typeof(max)}=min; pertubation::Vector=[])

Return the initial form of `f` with respect to valuation `val` and weight `w`.
If convention==min (default), it is computed in the min convention. If
convention==max, it is computed in the max convention.

For the definition of initial form in the min-convention, see
Section 2.4 of [MS15](@cite).

# Examples
```jldoctest
julia> Kxy, (x,y) = PolynomialRing(QQ,["x", "y"]);

julia> w = [1,1];

julia> val_2 = TropicalSemiringMap(QQ,2);

julia> val_trivial = TropicalSemiringMap(QQ);

julia> f = 2*x+2*y+1;

julia> initial(f,val_2,w)       # polynomial over GF(2)
1

julia> initial(f,val_trivial,w)
1
```
```jldoctest
julia> Kt,t = RationalFunctionField(QQ,"t");

julia> w = [1,1];

julia> Ktxy, (x,y) = PolynomialRing(Kt,["x", "y"]);

julia> f = t*x+t*y+1;

julia> val_t = TropicalSemiringMap(Kt,t);

julia> initial(f,val_t,w)       # polynomial over QQ
1
```
```jldoctest
julia> Kt,t = RationalFunctionField(GF(32003),"t");

julia> Ktxy, (x,y) = PolynomialRing(Kt,["x", "y"]);

julia> w = [1,1];

julia> f = t*x+t*y+1;

julia> val_t = TropicalSemiringMap(Kt,t);

julia> initial(f,val_t,w)       # polynomial over QQ
1
```
"""
function initial(f::MPolyElem, val::TropicalSemiringMap, w::Vector; pertubation::Vector=[])
  # compute the maximal weighted degrees
  # todo (optional):
  # currently, we iterate over the entire polynomial to compute the (terms with) maximal valuated weighted degrees
  # often this is not necessary as the polynomial is already sorted w.r.t. it
  if isempty(pertubation)
    vwd,vwds = valued_weighted_degree(f, val, w, return_vector=true)
  else
    vwd,vwdPerp,vwds,vwdsPerp = valued_weighted_degree(f, val, w, pertubation=pertubation, return_vector=true)
  end

  # initial(f) is the sum over all pi(c_alpha*t^-val(c_alpha))x^alpha
  # where c_alpha x^alpha is a term of maximal valued weighted degree
  # and pi is the map from the valued field to the residue field
  if is_valuation_trivial(val)
    t = val.valued_field(1)
  else
    t = val.valued_field(val.uniformizer_field)
  end
  kx, x = PolynomialRing(val.residue_field,[repr(x) for x in gens(parent(f))])
  R = val.valued_ring
  pi = val.residue_map

  initialf = MPolyBuildCtx(kx)
  if isempty(pertubation)
    for (vwdi,cf,expv) in zip(vwds,coefficients(f),exponent_vectors(f))
      if vwdi == vwd
        vcf = Int(val(cf),preserve_ordering=true)
        c = t^-vcf*cf   # make coefficient valuation 0
        cNum = numerator(c) # split up numerator and denominator as pi is only defined on the valued ring
        cDen = denominator(c)
        push_term!(initialf, pi(cNum)//pi(cDen), expv) # apply pi to both and divide the result
      end
    end
  else
    for (vwdi,vwdiPerp,cf,expv) in zip(vwds,vwdsPerp,coefficients(f),exponent_vectors(f))
      if vwdi == vwd && vwdiPerp == vwdPerp
        vcf = Int(val(cf),preserve_ordering=true)
        c = t^-vcf*cf
        cNum = numerator(c)
        cDen = denominator(c)
        push_term!(initialf, pi(cNum)//pi(cDen), expv)
      end
    end
  end

  return finish(initialf)
end
function initial(G::Vector, val::TropicalSemiringMap, w::Vector; pertubation::Vector=[])
  return [initial(g,val,w,pertubation=pertubation) for g in G]
end
export initial




@doc Markdown.doc"""
    initial(I::MPolyIdeal, val::TropicalSemiringMap, w::Vector; skip_groebner_basis_computation::Bool=false, skip_legality_check::Bool=false)

Return the initial ideal of `I` with respect to valuation `val` and weight `w`.
For the definition of initial ideal, see Section 2.4 of [MS15](@cite).

Use at your own risk: If `skip_groebner_basis_computation=true`, skips Groebner
basis computation. If `skip_legality_check=true`, skips check whether valuation
and weight vector are legal, i.e., if `I` is non-homogeneous, then `val` may
only be trivial and `w` may only have non-negative entries.

"""
function initial(I::MPolyIdeal, val::TropicalSemiringMap, w::Vector; skip_groebner_basis_computation::Bool=false)
  if !skip_groebner_basis_computation
    G = groebner_basis(I,val,w)
  else
    G = gens(G)
  end
  return ideal(initial(G,val,w))
end
