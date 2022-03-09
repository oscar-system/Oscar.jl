###
# Computing initial forms and initial ideals in Oscar
# ===================================================
###

@doc Markdown.doc"""
    valued_weighted_degree(f, val::ValuationMap, w::Vector; return_vector::Bool=false)

Returns the valued weighted degree of a polynomial `f` with respect to valuation `val` and weight vector `w`. In other words, returns the tropicalized polynomial of `f` with respect to valuation `val` evaluated at `w`.

If `return_vector=true`, returns a vector whose i-th entry is the valued weighted degree of the i-th term of `f`.

# Examples
```jldoctest
julia> Kxy, (x,y) = PolynomialRing(QQ,["x", "y"]);

julia> val_2 = ValuationMap(QQ,2);

julia> val_trivial = ValuationMap(QQ);

julia> w = [1,1];

julia> f = 2*x+2*y+1;

julia> valued_weighted_degree(f, val_2, w)

julia> valued_weighted_degree(f, val_2, w, return_vector=true)

julia> valued_weighted_degree(f, val_trivial, w, return_vector=true)

```
"""
function valued_weighted_degree(f::MPolyElem, val::ValuationMap, w::Vector; pertubation::Vector=[], return_vector::Bool=false)
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
  return vwd
end
export valued_weighted_degree



# # not wrong, but not sure whether needed
# function weighted_degree(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, w::Vector; return_vector::Bool=false)
#   trivial_val = ValuationMap(coefficient_ring(f))
#   return valued_weighted_degree(f, trivial_val, w, return_vector=return_vector)
# end
# export weighted_degree



@doc Markdown.doc"""
<<<<<<< HEAD
    initial(f, val::ValuationMap, w::Vector)
=======
    initial(f::MPolyElem, val::ValuationMap, w::Vector, convention::Union{typeof(min),typeof(max)}=min; pertubation::Vector=[])
>>>>>>> 16d4e1d18e (TropicalGeometry: renamed tropical_numbers to tropical_semiring)

Returns the initial form of `f` with respect to valuation `val` and weight `w`. If convention==min (default), it is computed in the min convention. If convention==max, it is computed in the max convention.

For the definition of initial form in the min-convention, see [Maclagan-Sturmfels, Section 2.4]

# Examples
```jldoctest
julia> Kxy, (x,y) = PolynomialRing(QQ,["x", "y"]);

julia> w = [1,1];

julia> val_2 = ValuationMap(QQ,2);

julia> val_trivial = ValuationMap(QQ);

julia> f = 2*x+2*y+1;

julia> initial(f,val_2,w)       # polynomial over GF(2)

julia> initial(f,val_trivial,w) # polynomial in the original ring

julia> initial(f,val_trivial,w,min) # same as above, as min-convention default

julia> initial(f,val_trivial,w,max) # different as above

julia> Kt,t = RationalFunctionField(QQ,"t");

julia> Ktxy, (x,y) = PolynomialRing(Kt,["x", "y"]);

julia> f = t*x+t*y+1;

julia> val_t = ValuationMap(Kt,t);

julia> initial(f,val_t,w)       # polynomial over QQ

julia> initial(f,val_trivial,w) # polynomial in the original ring

julia> Kt,t = RationalFunctionField(GF(32003),"t");

julia> Ktxy, (x,y) = PolynomialRing(Kt,["x", "y"]);

julia> f = t*x+t*y+1;

julia> val_t = ValuationMap(Kt,t)

julia> initial(f,val_t,w)       # polynomial over QQ

julia> initial(f,val_trivial,w) # polynomial in the original ring

```
"""
function initial(f, val::ValuationMap, w::Vector)
  # compute the maximal weighted degrees
  # todo (optional):
  # currently, we iterate over the entire polynomial to compute the (terms with) maximal valuated weighted degrees
  # often this is not necessary as the polynomial is already sorted w.r.t. it
  vwd,vwds = valued_weighted_degree(f, val, w, return_vector=true)

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
<<<<<<< HEAD
  for (vwdi,cf,expv) in zip(vwds,coefficients(f),exponent_vectors(f))
    if vwdi == vwd
      push_term!(initialf, pi(R(t^-val(cf)*cf)), expv)
=======
  if isempty(pertubation)
    for (vwdi,cf,expv) in zip(vwds,coefficients(f),exponent_vectors(f))
      if vwdi == vwd
        vcf = numerator(data(val(cf)))
        if convention(val)==max
          vcf = -vcf
        end
        c = t^-vcf*cf   # make coefficient valuation 0
        cNum = numerator(c) # split up numerator and denominator as pi is only defined on the valued ring
        cDen = denominator(c)
        push_term!(initialf, pi(cNum)//pi(cDen), expv) # apply pi to both and divide the result
      end
    end
  else
    for (vwdi,vwdiPerp,cf,expv) in zip(vwds,vwdsPerp,coefficients(f),exponent_vectors(f))
      if vwdi == vwd && vwdiPerp == vwdPerp
        vcf = numerator(data(val(cf)))
        if convention(val)==max
          vcf = -vcf
        end
        c = t^-vcf*cf
        cNum = numerator(c)
        cDen = denominator(c)
        push_term!(initialf, pi(cNum)//pi(cDen), expv)
      end
>>>>>>> 16d4e1d18e (TropicalGeometry: renamed tropical_numbers to tropical_semiring)
    end
  end

  return finish(initialf)
end
export initial


function initial(G::Vector{<:MPolyIdeal}, val::ValuationMap, w::Vector; skip_groebner_basis_computation::Bool=false)
  return [initial(g,val,w) for g in G]
end




@doc Markdown.doc"""
    initial(I::MPolyIdeal, val::ValuationMap, w::Vector; skip_groebner_basis_computation::Bool=false, skip_legality_check::Bool=false)

Returns the initial ideal of `I` with respect to valuation `val` and weight `w`. For the definition of initial ideal, see [Maclagan-Sturmfels, Section 2.4]

Use at your own risk: If `skip_groebner_basis_computation=true`, skips Groebner basis computation. If `skip_legality_check=true`, skips check whether valuation and weight vector are legal, i.e., if `I` is non-homogeneous, then `val` may only be trivial and `w` may only have non-negative entries.

# Examples
```jldoctest
julia> Kx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(QQ,6);

julia> Cyclic5Homogenized = ideal([x1+x2+x3+x4+x5,
                                   x1*x2+x2*x3+x3*x4+x1*x5+x4*x5,
                                   x1*x2*x3+x2*x3*x4+x1*x2*x5+x1*x4*x5+x3*x4*x5,
                                   x1*x2*x3*x4+x1*x2*x3*x5+x1*x2*x4*x5+x1*x3*x4*x5+x2*x3*x4*x5,
                                   -x0^5+x1*x2*x3*x4*x5]);

julia> Katsura5Homogenized = ideal([-x0+x1+2*x2+2*x3+2*x4+2*x5,
                                    -x0*x1+x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2,
                                    -x0*x2+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5,
                                    x2^2-x0*x3+2*x1*x3+2*x2*x4+2*x3*x5,
                                    2*x2*x3-x0*x4+2*x1*x4+2*x2*x5]);

julia> val_2 = ValuationMap(QQ,2); # 2-adic valuation

julia> val_3 = ValuationMap(QQ,3); # 3-adic valuation

julia> w = [0,0,0,0,0,0];

julia> initial(Cyclic5Homogenized, val_2, w)

julia> initial(Cyclic5Homogenized, val_3, w) # same as for val_2

julia> initial(Katsura5Homogenized, val_2, w)

julia> initial(Katsura5Homogenized, val_3, w) # different to val_2

julia> Kt,t = RationalFunctionField(QQ,"t");

julia> Ktx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(Kt,6);

julia> Cyclic5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Cyclic5Homogenized)]);

julia> Katsura5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Katsura5Homogenized)]);

julia> val_t = ValuationMap(Kt,t); # t-adic valuation

julia> initial(Cyclic5Homogenized_Kt, val_t, w) # same leading monomials as for val_2 and val_3

julia> initial(Katsura5Homogenized_Kt, val_t, w) # different leading monomials as for val_2
                                                 # same leading monomials as for val_3

```
"""
function initial(I::MPolyIdeal, val::ValuationMap, w::Vector; skip_groebner_basis_computation::Bool=false)
  if !skip_groebner_basis_computation
    G = groebner_basis(I,val,w)
  else
    G = gens(G)
  end
  return ideal(initial(G,val,w))
end
