###
# Computing initial forms and initial ideals in Oscar
# ===================================================
###

@doc Markdown.doc"""
    valued_weighted_degree()

Returns the valued weighted degree of a polynomial. The valued

# Examples
```jldoctest
julia> Kxy, (x,y) = QQ["x", "y"];

julia> val_2 = ValuationMap(QQ,2);

julia> f = 2*x+2*y+1;

julia> valued_weighted_degree(f,val_2,[1,1])
```
"""
function valued_weighted_degree(f, val, w; return_vector::Bool=false)
  # compute the weighted degrees
  vwds = [dot(w,alpha) for alpha in exponent_vectors(f)]
  # substract the coefficient valuations
  vwds -= [val(c) for c in coefficients(f)]

  # compute maximum
  vwd = max(vwds...)

  if return_vector
    return vwd,vwds
  end
  return vwd
end
export valued_weighted_degree



function weighted_degree(f::AbstractAlgebra.Generic.MPoly{<:RingElement}, w; return_vector::Bool=false)
  trivial_val = ValuationMap(coefficient_ring(f))
  return valued_weighted_degree(f, trivial_val, w, return_vector=return_vector)
end
export weighted_degree



@doc Markdown.doc"""
    initial(f,val,w)

Returns the initial form of `f` with respect to valuation `val` and weight `w`.

# Examples
```jldoctest
julia> Kxy, (x,y) = QQ["x", "y"];

julia> w = [1,1];

julia> val_2 = ValuationMap(QQ,2);

julia> val_trivial = ValuationMap(QQ);

julia> f = 2*x+2*y+1;

julia> initial(f,val_2,w)

julia> initial(f,val_trivial,w)

```
"""
function initial(f, val, w)
  # compute the maximal weighted degrees
  # todo (optional):
  # currently, we iterate over the entire polynomial to compute the (terms with) maximal valuated weighted degrees
  # often this is not necessary as the polynomial is already sorted w.r.t. it
  vwd,vwds = valued_weighted_degree(f, val, w, return_vector=true)

  # initial(f) is the sum over all pi(c_alpha*t^-val(c_alpha))x^alpha
  # where c_alpha x^alpha is a term of maximal valued weighted degree
  # and pi is the map from the valued field to the residue field
  if val.uniformizer==nothing
    t = val.valued_field(1)
  else
    t = val.valued_field(val.uniformizer)
  end
  kx, x = PolynomialRing(val.residue_field,[repr(x) for x in gens(parent(f))])
  pi = val.residue_map

  initialf = kx(0)
  for i in 1:length(vwds)
    if vwds[i] == vwd
      v = val(coeff(f,i))
      initialf += map_coefficients(pi,t^-v*term(f,i))
    end
  end

  return initialf
end
export initial



# function initial(I, val, w; skip_groebner_basis_computation::Bool=false, skip_legality_check::Bool=false)

#   if !skip_legality_check
#     check_legality(I,w,ignore_valuation=ignore_valuation,skip_groebner_basis_computation)
#   end

#   error("initial(ideal): current work in progress")
#   return I
# end

# returns true if the exponent vectors of g have the same sum
# return false otherwise
function sloppy_is_homogeneous(g)
  leadexpv,tailexpvs = Iterators.peel(exponent_vectors(g))
  d = sum(leadexpv)
  for tailexpv in tailexpvs
    if d != sum(tailexpv)
      return false
    end
  end
  return true
end


# todo: this needs a better name
# checks whether the following conditions are satisfied:
# - if weight vector has negative entries, then ideal needs to be homogeneous
# returns true if they are, returns false otherwise
function check_weight_vector_for_legality(I, val, w; skip_groebner_basis_computation::Bool=false)
  if !skip_groebner_basis_computation
    G = groebner_basis(I,complete_reduction=true)
  else
    G = gens(G)
  end

  is_ideal_homogeneous = true
  for g in G
    if sloppy_is_homogeneous(g)
      is_ideal_homogeneous = false
      break;
    end
  end

  K = coefficient_ring(G)
  if !is_ideal_homogeneous && is_valuation_trivial(val)
    error("check_legality: ideal needs to be homogeneous if computing w.r.t. non-trivial valuation")
  end

  is_weight_vector_nonnegative = true
  for wi in w
    if wi<0
      is_weight_vector_nonnegative = false
      break
    end
  end

  if !is_ideal_homogeneous && !is_weight_nonnegative
    error("check_legality: ideal needs to be homogenous if computing w.r.t. negative weight vector")
  end
end
