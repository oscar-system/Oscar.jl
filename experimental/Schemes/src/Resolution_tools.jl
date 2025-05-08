########################################################################
# intersection matrix of exceptional divisors -- surface case
########################################################################
@doc raw"""
    function intersection_matrix(phi::Union{BlowUpSequence,MixedBlowUpSequence})

Given a desingularization `phi` of a surface (2-dimensional reduced scheme), 
return the intersection matrix of the exceptional divisor of phi.

Return a tuple `M`, `v`, `M2`, `m` where 
`M` is the intersection matrix computed over the underlying field `k`, 
`v` is the list of k-irreducible components of the exceptional divisor ordered as in the columns of `M`, 
`M2` is the intersection matrix over the algebraic closure `kbar` of `k`,
`m` the number of kbar-irreducible components of each entry of `v`.

!!! warning
  This is only applicable, if the `phi` is the desingularization of a 2-dimensional reduced scheme, as the exceptional divisor of the desingularization of higher dimensional schemes is not a curve.

!!! note
  The intersection matrix referred to in textbooks is `M2`, as these usually restrict to the case of algebraically closed fields, but computations are usually performed over suitable subfields, e.g. `QQ` instead of `CC`. 

# Examples
```jldoctest
julia> R,(x,y,z) = polynomial_ring(QQ,3);

julia> W = AffineScheme(R);

julia> J = ideal(R,[x^2+y^2+z^5]);           ## A4 singularity, 2 pairs of abs. red. curves

julia> JS = IdealSheaf(W,J);

julia> Y = subscheme(JS);

julia> phi = desingularization(Y);

julia> L = intersection_matrix(phi);

julia> L[1]
[-4    2]
[ 2   -2]

julia> L[2]
2-element Vector{Oscar.AbsIdealSheaf}:
 Sheaf of prime ideals on scheme over QQ covered with 5 patches
 Sheaf of prime ideals on scheme over QQ covered with 5 patches

julia> L[3]
[-2    0    1    0]
[ 0   -2    0    1]
[ 1    0   -2    1]
[ 0    1    1   -2]

julia> L[4]
2-element Vector{Int64}:
 2
 2

```
"""
@attr Any function intersection_matrix(phi::Union{BlowUpSequence,MixedBlowUpSequence})
  phi.resolves_sing || error("intersection_matrix not available for partial desingularizations")
  !isdefined(phi, :is_embedded) || !phi.is_embedded || error("not available yet for embedded desingularization of curves")
  dim(domain(phi))==2 || error("not a surface -- exceptional locus not a graph")
  sl_orig = ideal_sheaf_of_singular_locus(codomain(phi))
  dim_sl_orig = dim(sl_orig)
  dim_sl_orig >= 0 || error("original scheme was non-singular, no exceptional curves added")
  dim_sl_orig == 0 || error("only available for isolated singularities")

## make sure that we have a strong resolution (i.e. exceptional divisor is simple normal crossing)
  phi = weak_to_strong_desingularization_surface(phi)

  cov_dom = simplified_covering(domain(phi))
  patches_scheme = patches(cov_dom)

## keep only non-empty exceptional curves
  ex_divs, dont_meet, caution_multi_charts = _cleanup_ex_div(phi)
## first determine intersection matrix: over given field
  inter_mat_k = zero_matrix(ZZ,length(ex_divs),length(ex_divs))

  # fill in the pairwise intersections
  # notice: by construction of the underlying resolution, it suffices to look at one chart
  #         in which a non-empty intersection of exceptional k-components manifests itself
  for i in 1:nrows(inter_mat_k)
    for j in i+1:nrows(inter_mat_k)
      inter_id = one(OO(patches_scheme[1]))
      !((i,j) in dont_meet) ||  continue                # cannot meet, entry stays 0
      if !((i,j) in caution_multi_charts)
        found_index = findfirst(V->(!is_one(ex_divs[i](V) + ex_divs[j](V))), patches_scheme)
        found_index !== nothing || continue
        U = patches_scheme[found_index]
        inter_id = ex_divs[i](U) + ex_divs[j](U)
        inter_mat_k[i,j] = vector_space_dimension(quo(base_ring(inter_id),inter_id)[1])
      else
        temp_inter = ex_divs[i] + ex_divs[j]
        !is_one(temp_inter) || continue
        tempint = 0
        for U in patches_scheme
          inter_id = temp_inter(U) + decomposition_info(U)
          tempint += vector_space_dimension(quo(base_ring(inter_id),inter_id)[1])
        end
        inter_mat_k[i,j] = tempint
      end
    end
  end

  # and fill in the lower triangular part as the property 'intersects' is symmetric
  inter_mat_k = inter_mat_k + transpose(inter_mat_k)

## get ready to compute self intersection numbers:
## choose curve passing through component of singular locus and decompose its full preimage
## (we need the strict transform and the exceptional multiplicities... --> ex_mult_dict )
  decomp_sl_orig = maximal_associated_points(sl_orig)   # we might have several singular points

  # we need a curve passing through a component of the singular locus, to find the self intersection numbers
  # of the exceptional curves arising from this component of the singular locus
  patches_orig = patches(default_covering(scheme(sl_orig)))
  ex_mult_dict = IdDict{AbsIdealSheaf,Tuple{AbsIdealSheaf,Vector{Int}}}()
  for I in decomp_sl_orig
    found_index = findfirst(V -> !is_one(I(V)), patches_orig)
    U = patches_orig[found_index]
    dim(I) == 0 || error("case of non-isolated singularities not implemented yet")
    if vector_space_dimension(quo(OO(U),I(U))[1]) == 1
      a = rational_point_coordinates(saturated_ideal(I(U)))
      l = findfirst(x -> is_prime(ideal(OO(U),gen(OO(U),x)-a[x])),1:ngens(OO(U)))
      if l !== nothing
        p = gen(OO(U),l) - a[l]
      else
        p = sum(gen(OO(U),x) - a[x] for x in 1:ngens(OO(U)))
        if !is_prime(ideal(OO(U),[p]))
          decomp_h = minimal_primes(ideal(OO(U),[p]))
          h_ind = findfirst(x -> is_subset(x,I(U)), decomp_h)
          h_ind !== nothing || error("choose hypersurface through point: no irreducible one found")
          p = small_generating_set(decomp_h[h_ind])
        end
      end
    else
      l = findfirst(x -> (deg(x) == 1 && is_prime(ideal(OO(U),[x]))), gens(I(U)))
      l !== nothing || error("choose hypersurface: random choice not implemented yet")
      p = gen(I(U),l)
    end
    H = IdealSheaf(scheme(sl_orig),U,ideal(OO(U),p))
    full_preimage_H = pullback(phi,H)
    H_strict = strict_transform(phi,H)
    decomp_H = _decompose_pullback_simple(full_preimage_H,ex_divs)
    ex_mult_dict[H] = (H_strict,decomp_H)
  end

## use ex_mult_dict to determine self intersection numbers from the fact that intersection numbers
## in codomain(phi) and the corresponding ones of the full preimages in domain(phi) are equal
  for (H,temp) in ex_mult_dict
    H_strict = temp[1]
    v = temp[2]
    # determine first those E_i, whose self intersection number we can determined with this H,
    # ignoring the ones which have already been treated
    # (Background; different singular components lead to disjoint parts of the exceptional locus,
    #              not every H necessarily meets all of these parts)
    E_todo = findall(i -> (inter_mat_k[i,i] == 0 && v[i] != 0), 1:length(v))
    for i in E_todo
      strict_inter = H_strict + ex_divs[i]
      # put the length of the zero-dimensional scheme strict_inter, i.e. H_strict . E_i
      #  into strict_summand
      found_index = findfirst(V -> !is_one(strict_inter(V)), patches_scheme)
      strict_summand = (found_index == nothing ? 0 :
           vector_space_dimension(quo(base_ring(OO(patches_scheme[found_index])),
                 saturated_ideal(strict_inter(patches_scheme[found_index]))
                       + modulus(OO(patches_scheme[found_index])))[1]))
      # next we consider H . equidimensional_hull(sl_orig)  -- zero, if dim(sl_orig) == 0
      dim_sl_orig == 0 || error("not implemented yet")
      orig_inter = 0
      divtemp = divides(orig_inter
                       - sum([m*e for (m,e) in zip(v,[inter_mat_k[i,l] for l in 1:nrows(inter_mat_k)])])
                       - strict_summand,
                     v[i])
      divtemp[1] == true || error("Cannot happen by theory: not divisible!")
      inter_mat_k[i,i] = divtemp[2]
    end
  end

  restemp = _pass_to_kbar_intermat(inter_mat_k, ex_divs)
  result = inter_mat_k, ex_divs, restemp[1], restemp[2]
  return result
end

function _pass_to_kbar_intermat(M::ZZMatrix, ex_divs::Vector{AbsIdealSheaf})
  ncols(M) == nrows(M) || error("not a square matrix")
  ncols(M) == length(ex_divs) || error("divisor list length does not match matrix size")

## these components are certainly absolutely irreducible: intersection multiplicity 1
  abs_irred_list = filter(a -> any(b -> (M[a,b] == 1 || M[b,a] ==1), 1:ncols(M)),1:ncols(M))

## now the reducible ones -- infer as much as possible to avoid absolute_primary_decomposition
  abs_data = Vector{Tuple{Int,Int,Int}}()  # (age of divisor,#abs_components,#A1-points)
  for a in 1:ncols(M)
    if a in abs_irred_list
    ## irreducible ==> all data known directly
      push!(abs_data,(a,1,0))
      continue
    end
    ## could still be absolutely irreducible, but also absolutely reducible
    ## first check, whether we can deduce that it is absolutely reducible
    b = findfirst(c -> (M[c,a] > 1 || M[a,c] > 1), abs_irred_list)
    if b !== nothing
      ## inferable from other data
      if has_attribute(ex_divs[a],:A1count)
        ## we already counted A1s
        c = get_attribute(ex_divs[a],:A1count)
        push!(abs_red_data, (a,M[b,a],c))
      else
        ## we still need to count A1s
        worse_sing,A1count = curve_sing_A1_or_beyond(ex_divs[a])
        is_one(worse_sing) || error("original desingularization was not strong desingularization")
        push!(abs_data, (a,M[b,a],A1count))
      end
      continue
    end
    ## we ended up in the default: not inferable, fill in A1s, but postpone absolute primary decomp.
    worse_sing,A1count = curve_sing_A1_or_beyond(ex_divs[a])
    is_one(worse_sing) || error("original desingularization was not strong desingularization")
    push!(abs_data, (a,-1,A1count))   ## negative number to be corrected later
  end

  ## now a second round of trying to infer instead of absolute_primary_decomposition
  patches_scheme = patches(simplified_covering(scheme(ex_divs[1])))
  todo_list = findall(a -> a[2] < 0, abs_data)
  for i in todo_list
    ## try intersections with known data
    j = findfirst(a->((M[a,i]>0 || M[i,a]>0) && abs_data[a][2] >0),1:ncols(M))
    if (j !== nothing) && (M[i,j] > abs_data[j][2])
      ## either the k_bar components meet pairwise or all meet all
      ## ==> ex_divs[i] not absolutely irreducible and all meet all
        is_sane, mult_i = divides(M[i,j],abs_data[j][2])
        is_sane || error("inconsistency of multiplicities of divisors over k_bar")
        abs_data[i][2] = mult_i
        continue
    end
    ## all else failed, we cannot avoid absolute primary decomposition
    found_index = findfirst(V -> !is_one(ex_divs[i](V)), patches_scheme)
    U = patches_scheme[found_index]
    decomp = absolute_primary_decomposition(ex_divs[i](U))
    length(decomp) == 1 || error("problem with decomposition")
    abs_data[i] = (abs_data[i][1],decomp[1][4],abs_data[i][3])
  end

  ## don't go any further, if we are done
  findfirst(v -> v[2]!=1 , abs_data) !== nothing || return M,[1 for i in 1:nrows(M)]

  ## fill the non-diagonal entries of the intersection matrix over algebraic closure
  ncols_kbar = sum([a[2] for a in abs_data])
  inter_mat = zero_matrix(ZZ,ncols_kbar,ncols_kbar)

  h_position = 0                                     # start in first row of M
  for i in 1:length(abs_data)
    v_position = h_position + abs_data[i][2]         # start after current block
    ## iterate over rows of original matrix leading to blocks of inter_mat
    for j in i+1:length(abs_data)
      ## iterate over the columns right of the diagonal
      ## M[i,j] == 0 only contributes zero matrix and can thus be omitted
      if M[i,j] == 0
        ## do nothing -- just for easier debugging
      elseif abs_data[i][2]*abs_data[j][2] == M[i,j]
        ## all meet all
        inter_mat[(h_position+1):(h_position+abs_data[i][2]),(v_position+1):(v_position+abs_data[j][2])] =
           ZZMatrix(ones(Int,abs_data[i][2],abs_data[j][2]))
      elseif abs_data[i][2] == M[i,j] && abs_data[j][2] == M[i,j]
        ## they meet pairwise
        inter_mat[h_position+1:h_position+abs_data[i][2],v_position+1:v_position+abs_data[j][2]] =
           identity_matrix(ZZ,abs_data[i][2])
      else
        ## should not happen -- error message for tracking down problems
        error("mismatch of divisor multiplicities and intersection multiplicity")
      end
      v_position += abs_data[j][2]      # continue one block to the right
    end
    h_position += abs_data[i][2]        # continue with subsequent row
  end

  inter_mat = inter_mat + transpose(inter_mat)

  ## now the self intersection numbers
  h_position = 0
  for i in 1:length(abs_data)
    is_sane, self_intersection = divides(- 2*abs_data[i][3] + M[i,i], abs_data[i][2])
    is_sane || error("problem: self intersection kbar went wrong")
    ## first summand self-intersection
    inter_mat[(h_position+1):(h_position+abs_data[i][2]),(h_position+1):(h_position+abs_data[i][2])] +=
           self_intersection * identity_matrix(ZZ,abs_data[i][2])
    ## second summand pairwise intersection of k_bar irreducible components
    if abs_data[i][3] != 0
      ## either they do not meet at all or all meet all as they, so here all meet all
      ## hence we put a 1 everywhere in the block, but not on the diagonal of the block
      inter_mat[(h_position+1):(h_position+abs_data[i][2]),(h_position+1):(h_position+abs_data[i][2])] +=
        ZZMatrix(ones(Int,abs_data[i][2],abs_data[i][2])) - identity_matrix(ZZ,abs_data[i][2])
    end
    h_position += abs_data[i][2]
  end

  return inter_mat,[a[2] for a in abs_data]

end

function _cleanup_ex_div(phi::BlowUpSequence)
  return _cleanup_ex_div(mixed_blow_up_sequence(phi))
end

function _cleanup_ex_div(phi::MixedBlowUpSequence)
## initialization
  ex_divs = phi.ex_div
  ret_divs = Vector{AbsIdealSheaf}()
  dont_meet_raw = ( isdefined(phi, :dont_meet) ? phi.dont_meet : Tuple{Int,Int}[])
  dont_meet = Vector{Tuple{Int,Int}}()
  caution_multi_charts_raw = ( isdefined(phi, :caution_multi_charts) ? phi.caution_multi_charts : Tuple{Int,Int}[])
  caution_multi_charts = Vector{Tuple{Int,Int}}()
  skip_list = Vector{Int}()
  skip_count = length(skip_list)
  offset = 0
  remember_orig = Vector{Int}()

## throw away components which have disappeared due to blow-ups
  div_list_raw = Vector{AbsIdealSheaf}()
  for i in 1:length(ex_divs)
    if dim(ex_divs[i]) == 1
      # divisor still relevant
      push!(div_list_raw,ex_divs[i])
    else
      # divisor obsoleted by later blow-up
      push!(skip_list,i)
      skip_count = skip_count+1
    end
  end

  # correct dont_meet and caution_multi_charts
  if skip_count > 0
    for i in 0:skip_count-1
      dont_meet_raw = [(a >= skip_list[skip_count-i] ? a-1 : a,
                        b >= skip_list[skip_count-i] ? b-1 : b) for (a,b) in dont_meet_raw]
      dont_meet = [a for a in dont_meet_raw if a[1] > 0]
      caution_multi_charts_raw = [(a >= skip_list[skip_count-i] ? a-1 : a,
                        b >= skip_list[skip_count-i] ? b-1 : b)
                        for (a,b) in caution_multi_charts_raw]
      caution_multi_charts = [a for a in caution_multi_charts_raw if a[1] > 0]
    end
  end

  # decompose divisors over given base field
  for i in 1:length(div_list_raw)
    I_list = maximal_associated_points(div_list_raw[i])
    append!(remember_orig, [i for j in 1:length(I_list)])   # note the original index
    append!(ret_divs,I_list)
  end

  if length(ret_divs) > length(div_list_raw)
    for (i,j) in dont_meet_raw
      append!(dont_meet,[(a,b) for a in findall(==(i),remember_orig) for b in findall(d -> d == j, remember_orig)])
      append!(caution_multi_charts,[(a,b) for a in findall(==(i),remember_orig) for b in findall(d -> d == j, remember_orig)])
    end
  end

  return ret_divs, dont_meet, caution_multi_charts

end

function _decompose_pullback_simple(H::AbsIdealSheaf, ex_divs::Vector{AbsIdealSheaf})
  @assert scheme(ex_divs[1]) == scheme(H)
  patches_H = patches(simplified_covering(scheme(H)))
  v = Vector{Int}()
  for E in ex_divs
    found_index = findfirst(V -> (!is_one(E(V)) && !is_one(H(V))), patches_H)
    if found_index == nothing
      push!(v,ZZ(0))
      continue
    end
    U = patches_H[found_index]
    _, v_cur = iterated_quotients(H(U), E(U),0)
    push!(v,v_cur)
  end
  return v
end
