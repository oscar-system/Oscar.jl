###
# Computing (tropical) Groebner bases in Oscar
# ============================================
#
# For a definition of tropical Groebner basis see Section 2.4 in:
#   D. Maclagan, B. Sturmfels: Introduction to tropical geometry
# To see how they can be computed using standard bases see:
#   T. Markwig, Y. Ren: Computing tropical varieties over fields with valuation
###


#=======
tropical Groebner basis
Example:
val_2 = ValuationMap(QQ,2)
Kx,(x,y,z) = PolynomialRing(QQ,3)
w = [0,0,0]
I = ideal([x+2*y,y+2*z])
groebner_basis(I,val_2,w)

Kt,t = RationalFunctionField(QQ,"t")
val_t = ValuationMap(Kt,t)
Ktx,(x,y,z) = PolynomialRing(Kt,3)
w = [0,0,0]
I = ideal([x+t*y,y+t*z])
groebner_basis(I,val_t,w,return_initial=true)
=======#
@doc Markdown.doc"""
    groebner_basis(I::Ideal, val::ValuationMap, w::Vector)

Computes a Groebner basis of `I` over a field with valuation `val` with respect to weight vector `w`, that is a finite generating set of `I` whose initial forms generate the initial ideal with respect to `w`.

For the definitions of initial form, initial ideal and Groebner basis see [Maclagan-Sturmfels, Section 2.4].

# Warning
`I` must be homogeneous if `val` is non-trivial or `w` contains non-positive entries. If `val` is trivla and `w` contains only non-negative entries, then what is computed is a regular Groebner basis with respect to a weighted ordering with weight vector `w`.

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

julia> groebner_basis(Cyclic5Homogenized, val_2, w)

julia> groebner_basis(Cyclic5Homogenized, val_3, w) # same as for val_2

julia> groebner_basis(Katsura5Homogenized, val_2, w)

julia> groebner_basis(Katsura5Homogenized, val_3, w) # different to val_2

julia> Kt,t = RationalFunctionField(QQ,"t");

julia> Ktx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(Kt,6);

julia> Cyclic5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Cyclic5Homogenized)]);

julia> Katsura5Homogenized_Kt = ideal([change_coefficient_ring(Kt,f) for f in gens(Katsura5Homogenized)]);

julia> val_t = ValuationMap(Kt,t); # t-adic valuation

julia> groebner_basis(Cyclic5Homogenized_Kt, val_t, w) # same leading monomials as for val_2 and val_3

julia> groebner_basis(Katsura5Homogenized_Kt, val_t, w) # different leading monomials as for val_2
                                                        # same leading monomials as for val_3
```
"""
function groebner_basis(I::MPolyIdeal,val::ValuationMap,w::Vector{<: Union{Int,Rational{Int},fmpz,fmpq} }; pertubation::Vector=[], skip_legality_check::Bool=false)

  ###
  # Step 0: check legality of input unless stated otherwise
  #   If val is non-trivial or w is not non-negative, I must be homogeneous
  #   This is because otherwise the localization with respect to the non-global ordering will change the ideal
  ###
  if !skip_legality_check
    check_legality(I, val, w)
  end

  ###
  # Step 1: Compute a standard basis in the simulation ring
  ###
  vvI = simulate_valuation(I,val)
  Rtx = base_ring(vvI)
  if isempty(pertubation)
    w = simulate_valuation(w,val)
    S,_ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(w)*Singular.ordering_dp())
  else
    w,u = simulate_valuation(w,pertubation,val)
    S,_ = Singular.PolynomialRing(singular_ring(base_ring(Rtx)), map(string, Nemo.symbols(Rtx)), ordering = Singular.ordering_a(w)*Singular.ordering_a(u)*Singular.ordering_dp())
  end
  SI = Singular.Ideal(S, [S(g) for g in gens(vvI)])
  vvGB = Singular.gens(Singular.satstd(SI,Singular.MaximalIdeal(S,1)))

  ###
  # Step 2: tighten simulation and return
  ###
  vvGB = [tighten_simulation(Rtx(g),val) for g in vvGB]
  return [g for g in desimulate_valuation(vvGB,val) if !iszero(g)]

end


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


# checks whether the ideal is homogeneous if val is non-trivial or w has negative entries
function check_legality(I::MPolyIdeal, val::ValuationMap, w::Vector)

  is_weight_vector_nonnegative = true
  for wi in w
    if wi<0
      is_weight_vector_nonnegative = false
      break
    end
  end

  if is_valuation_nontrivial(val) || is_weight_vector_nonnegative
    for g in gens(I) # todo: interreduce generators before test for homogeneity
      if !sloppy_is_homogeneous(g)
        error("ideal needs to be homogeneous if computing w.r.t. non-trivial valuation")
      end
    end
  end

end


###
# Returns a reduced GB if given a GB.
# Note that a reduced GB is not always desirable as its simulation might have a smaller leading ideal
# (due to possibly higher powers of t)
# This is why reducing a GB is generally not recommended, unless it is absolutely necessary.
# And the reduced GB should be discarded for the unredued GB as soon as it is of no use.
###
#=======
tropical Groebner basis reduction
Example:
Kx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(QQ,6);
Cyclic5Homogenized = ideal([x1+x2+x3+x4+x5,
                            x1*x2+x2*x3+x3*x4+x1*x5+x4*x5,
                            x1*x2*x3+x2*x3*x4+x1*x2*x5+x1*x4*x5+x3*x4*x5,
                            x1*x2*x3*x4+x1*x2*x3*x5+x1*x2*x4*x5+x1*x3*x4*x5+x2*x3*x4*x5,
                            -x0^5+x1*x2*x3*x4*x5]);
Katsura5Homogenized = ideal([-x0+x1+2*x2+2*x3+2*x4+2*x5,
                             -x0*x1+x1^2+2*x2^2+2*x3^2+2*x4^2+2*x5^2,
                             -x0*x2+2*x1*x2+2*x2*x3+2*x3*x4+2*x4*x5,
                             x2^2-x0*x3+2*x1*x3+2*x2*x4+2*x3*x5,
                             2*x2*x3-x0*x4+2*x1*x4+2*x2*x5]);
val = ValuationMap(QQ,2);
w = [0,0,0,0,0,0];
G = groebner_basis(Cyclic5Homogenized, val, w)
# G = groebner_basis(Katsura5Homogenized, val, w)
interreduce_tropically(G,val,w)


Ks,s = RationalFunctionField(QQ,"s");
Ksx,(x0,x1,x2,x3,x4,x5) = PolynomialRing(Ks,6);
Cyclic5Homogenized_Ks = ideal([change_coefficient_ring(Ks,f) for f in gens(Cyclic5Homogenized)]);
Katsura5Homogenized_Ks = ideal([change_coefficient_ring(Ks,f) for f in gens(Katsura5Homogenized)]);
val = ValuationMap(Ks,s); # t-adic valuation
w = [0,0,0,0,0,0];
G = groebner_basis(Cyclic5Homogenized_Ks, val, w)
# G = groebner_basis(Katsura5Homogenized_Ks, val, w)

=======#
function interreduce_tropically(G::Vector{<:MPolyElem}, val::ValuationMap, w::Vector; pertubation::Vector=[])

  # println("================================================== inside interreduce_tropically")
  # println(G)
  # println(w)
  # println(pertubation)

  ###
  # Step 0: simulate valuation and change coefficient ring to valued field
  ###
  vG = simulate_valuation(G,val,coefficient_field=true)
  Rtx = parent(vG[1])
  if isempty(pertubation)
    vw = simulate_valuation(w,val)
    S,_ = Singular.PolynomialRing(singular_ring(val.valued_field),
                                  map(string, Nemo.symbols(Rtx)),
                                  ordering = Singular.ordering_a(vw)*Singular.ordering_dp())
  else
    vw,vu = simulate_valuation(w,pertubation,val)
    S,_ = Singular.PolynomialRing(singular_ring(val.valued_field),
                                  map(string, Nemo.symbols(Rtx)),
                                  ordering = Singular.ordering_a(vw)*Singular.ordering_a(vu)*Singular.ordering_dp())
  end
  sG = [S(change_base_ring(val.valued_field,g)) for g in vG] # todo: remove workaround when fixed


  ###
  # Step 1: seperate elements of sG by degree in x
  ###
  sG_degrees = [x_degree(sg) for sg in sG]
  sG_slices = [[] for d in 1:max(sG_degrees...)]
  sG_slice0 = []
  for (sg,d) in zip(sG,sG_degrees)
    if d>0
      push!(sG_slices[d],sg)
    else
      push!(sG_slice0,sg)
    end
  end

  if length(sG_slice0)!=1
    error("input simulated Groebner basis suspicious number of x-degree 0 elements")
  end


  ###
  # Step 2: sort and interreduce_tropically each degree slice
  ###
  Singular.libSingular.set_option("OPT_INFREDTAIL", true)
  for (d,H) in enumerate(sG_slices)

    # skip H if it is empty
    if (isempty(H))
      continue
    end

    # println("================================================== H start of step 2")
    # display(H)

    # sort H by leading monomial from large to small
    sort!(H,rev=true)

    # println("================================================== H after sort")
    # display(H)

    # first pass, remove leading x-monomial of H[i] from H[j] for i<j
    for i in 1:length(H)-1
      for j in i+1:length(H)
          H[j] = Singular.reduce(H[j],Singular.std(Singular.Ideal(S,H[i])))
      end
    end

    # println("================================================== H after 1st pass")
    # display(H)

    # second pass, remove leading x-monomial of H[j] from H[i] for i<j
    for i in 1:length(H)-1
      for j in i+1:length(H)
          H[i] = Singular.reduce(H[i],Singular.std(Singular.Ideal(S,H[j])))
      end
    end

    # println("================================================== H after 2nd pass")
    # display(H)

    # overwrite old slice with the new reduced slice
    sG_slices[d] = H
  end


  ###
  # Step 3: reduce each degree slice by its lower-degree predecessors
  ###
  # records whether we have encountered the lowest degree slice
  # this is so we can skip the lowest degree slice without any reduction
  sG_reduced = []
  for H in sG_slices

    # if slice is empty, skip to next slice
    if isempty(H)
      continue
    end

    # if no reducers exist, add current slice to reducers and skip to next slice
    if isempty(sG_reduced)
      append!(sG_reduced,H)
      continue
    end

    # H will increase in size during the reduction process
    # the following vector keeps track which of its entries are original elements
    H_original_poly_indices = collect(1:length(H))
    # records up to which term each original element of H has been reduced
    # * 1 means only leading monomial is reduced which is always true
    # * 0 means that all elements of H have been reduced
    H_original_tail_indices = ones(Int,length(H))
    zero_vector = zeros(Int,length(H))
    while H_original_tail_indices>zero_vector
      # Step 3.1: Find the maximal monomial to be reduced
      monomials_to_be_reduced = [first(Iterators.drop(monomials(H[i]),j))
                                 for (i,j) in zip(H_original_poly_indices,
                                                  H_original_tail_indices) if j>0]
      index_to_be_reduced = argmax(monomials_to_be_reduced)
      monomial_to_be_reduced = monomials_to_be_reduced[index_to_be_reduced]

      # Step 3.2: Find the reducer that is maximal once multiplied to the right x-degree
      #   as the x-monomial will be the same as that of monomial_to_be_reduced,
      #   this is equivalent to finding the reducer with maximal t-degree.
      #   (the maximality is optional, but comes with several practical advantages)
      new_reducer_t_degree = -1
      new_reducer = 0
      for g in sG_reduced
        b,q = divides(monomial_to_be_reduced,leading_monomial(g))
        if b && new_reducer_t_degree<Singular.leading_exponent_vector(g)[1]
          new_reducer = q*g
          new_reducer_t_degree = Singular.leading_exponent_vector(g)[1]
        end
      end

      # if there is no reducer, increment H_original_tail_indices_to_be_reduced
      #   and skip to next monomial_to_be_reduced
      if iszero(new_reducer)
        H_original_tail_indices[index_to_be_reduced] += 1
        # set H_original_tail_indices to 0 if length exceeded
        if length(H[H_original_poly_indices[index_to_be_reduced]])>=H_original_tail_indices[index_to_be_reduced]
          H_original_tail_indices[index_to_be_reduced] = 0
        end
        continue
      end

      # Step 3.3: Add new_reducer to H and update H_original_element_indices
      new_reducer_position = searchsortedfirst(H,new_reducer)
      insert!(H,new_reducer_position,new_reducer)
      for i in 1:length(H_original_poly_indices)
        if H_original_poly_indices[i] >= new_reducer_position
          H_original_poly_indices[i] += 1
        end
      end

      # Step 3.4: Reduce new H
      # first pass, remove leading x-monomial of H[i] from H[j] for i<j when new_reducer_position<=j
      for i in 1:length(H)-1
        for j in max(new_reducer_position,i+1):length(H)
          H[j] = Singular.reduce(H[j],Singular.std(Singular.Ideal(S,H[i])))
        end
      end
      # second pass, remove leading x-monomial of H[j] from H[i] for i<j when new_reducer_position<=j
      for i in 1:length(H)-1
        for j in max(new_reducer_position,i+1):length(H)
          H[i] = Singular.reduce(H[i],Singular.std(Singular.Ideal(S,H[j])))
        end
      end

      # Step 3.5: update H_original_tail_indices by replacing larger equal length entries with 0
      for (i,j) in enumerate(H_original_poly_indices)
        if H_original_tail_indices[i]>=length(H[j])
          H_original_tail_indices[i] = 0
        end
      end
    end

    # append original elements to sG_reduced
    append!(sG_reduced,H[H_original_poly_indices])
  end
  Singular.libSingular.set_option("OPT_INFREDTAIL", false)


  ###
  # Step 4: return reduced GB
  ###
  vG = [Rtx(sg) for sg in sG_reduced]
  return desimulate_valuation(vG,val)

  ###
  # old code:
  ###
  # if is_valuation_trivial(val)
  #   # todo: just call Singular.interred using the correct ordering
  # else
  #   sort!(vvGB,lt=x_monomial_lt) # sort vvGB by their leading x monomial from small to large
  #   Singular.libSingular.set_option("OPT_INFREDTAIL", true)
  #   for i in 1:length(vvGB)-1
  #     for j in i+1:length(vvGB)
  #       t_ecart = x_monomial_ecart(vvGB[j],vvGB[i])
  #       if t_ecart>=0
  #         vvGB[j] = Singular.reduce(val.uniformizer_ring^t_ecart*vvGB[j],Singular.std(Singular.Ideal(S,vvGB[i])))
  #         vvGB[j] = S(tighten_simulation(Rtx(vvGB[j]),val))
  #       end
  #     end
  #   end
  #   Singular.libSingular.set_option("OPT_INFREDTAIL", false)
  # end

end
export interreduce_tropically

###
# returns true if (leading x-monomial of f) <_lex (leading x-monomial of g)
# returns false otherwise
###
function x_degree(f::Singular.spoly)
  exp_f = Singular.leading_exponent_vector(f)
  return sum(exp_f)-exp_f[1]
end
export x_degree

# ###
# # returns true if (leading x-monomial of f) <_lex (leading x-monomial of g)
# # returns false otherwise
# ###
# function x_monomial_lt(f::Singular.spoly, g::Singular.spoly)
#   exp_x_f = copy(Singular.leading_exponent_vector(f))
#   exp_x_g = copy(Singular.leading_exponent_vector(g))
#   popfirst!(exp_x_f)
#   popfirst!(exp_x_g)
#   return exp_x_f<exp_x_g
# end
# export x_monomial_lt


# ###
# # returns true if x^expv_g divides x^expv_f
# # returns false otherwise
# ###
# function x_monomial_divides(exp_x_g::Vector,exp_x_f::Vector)
#   for (eg,ef) in zip(exp_x_g,exp_x_f)
#     if eg>ef
#       return false
#     end
#   end
#   return true
# end
# export x_monomial_divides


# ###
# # if the leading x-monomial of g divides x-monomials in f
# #   returns l=max(0, (t_g-exponent of g) - (t-exponents of f))
# #   so that f*t^l can be reduced by g to eliminate all the x-monomials
# # otherwise, returns -1
# ###
# function x_monomial_ecart(f::Singular.spoly, g::Singular.spoly)
#   exp_x_g = copy(Singular.leading_exponent_vector(g))
#   exp_t_g = popfirst!(exp_x_g)
#   e = 0
#   dividend_found = false
#   for exp_f in exponent_vectors(f)
#     exp_x_f = copy(exp_f)
#     exp_t_f = popfirst!(exp_x_f)
#     if x_monomial_divides(exp_x_g,exp_x_f)
#       e = max(e,exp_t_g-exp_t_f)
#       dividend_found = true
#     end
#   end
#   if dividend_found
#     return e
#   end
#   return -1
# end
# export x_monomial_ecart
