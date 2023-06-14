function has_du_val_singularities(F::MPolyDecRingElem)::Bool

  R = parent(F)
  is_graded(R) || error("only available for graded rings")
  is_standard_graded(R) || error("only available for standard grading")
  ncharts = ngens(R)

  ncharts == 4 || error("not a surface")            ## hypersurface in P^3

  JF = jacobi_ideal(F)
  I_sl = JF + ideal(R,F)
  dim(I_sl) > 0 || return true                      ## smooth, hence 'at most duVal'
  decomp = minimal_primes(I_sl)

  not_over_k = Vector{Tuple{MPolyIdeal,Int}}()
  dehom_dict = IdDict{Int,Ideal}()

  for I in decomp
    d = dim(I)
    d < 2 || return false                           ## non-isolated is not duVal
    
    ## It is definitely a point now, i.e. I becomes maximal after dehomog.
    ## It shows up in at least one dehomogenizatied chart
    i = 0
    I_dehomog = dehomogenization(I,1)               ## initialize I_dehomog
    while i < ncharts
      i=i+1
      I_dehomog = dehomogenization(I,i)
      is_one(I_dehomog) || break
    end

    ## cache dehomogenized ideal
    if !haskey(dehom_dict,i)
      dehom_dict[i] = dehomogenization(I_sl,i)
    end

    ## make sure to do only the geometrically irreducible points first
    ## and postpone all geometrically reducible points
    if vdim(quo(base_ring(I_dehomog),I_dehomog)[1]) > 1 
      push!(not_over_k, (I_dehomog,i))
      continue
    end

    ## now do the check for this geometrically irreducible point
    _check_duval_at_point(I_dehomog,dehom_dict[i]) || return false    

  end

  ## handle geometrically reducible points afterwards
  for (I,i) in not_over_k
    li = absolute_primary_decomposition(I)
    for (_, _, I2, _) in li
      ## prepare setting to check over algebraic extension
      r_changed = base_ring(I2)
      kk = coefficient_ring(r_changed)
      I_sl_change = ideal(r_changed, [change_coefficient_ring(kk,a) for a in gens(dehom_dict[i])])
      ## now do the check (at one of the conjugate points)
      _check_duval_at_point(I2,I_sl_change) || return false
    end
  end

  return true
end

function _check_duval_at_point(I_r::MPolyIdeal,I_sl_r::MPolyIdeal)::Bool
    ## initialize
    r = base_ring(I_r)

    ## localize at point
    a = Oscar.rational_point_coordinates(I_r)
    U = complement_of_point_ideal(r,a)
    rl, loc_map = localization(r,U)
    I_loc = loc_map(I_sl_r)

    o = negdegrevlex(gens(r))

    ## do quick version of beginning of Arnold's classifier
    tau=vdim(quo(rl,I_loc)[1])

    corank = vdim(quo(rl,I_loc+loc_map(I_r)^2)[1]) - 1
                                         # only count degree 1 monomials
    corank < 3 || return false           # at least T_3,3,3, not duVal
    corank > 1 || return true            # A series

    cubiccount = vdim(quo(rl,I_loc+loc_map(I_r)^3)[1]) - corank - 1
                                         # only count degree 2 monomials
    cubiccount < 3 || return false       # at least X_9, not duVal 
    cubiccount > 1 || return true        # D series

    tau < 9 || return false              # at least J_10           
    return true                          # E_k, k \in {6,7,8}
end

function has_du_val_singularities(J::MPolyIdeal{<:MPolyDecRingElem})
  R = base_ring(J)
  is_graded(R) || error("only available for graded rings")
  is_standard_graded(R) || error("only available for standard grading")
  ncharts = ngens(R)

  # mind that dimensions are shifted by one as the methods refer to
  # the polynomial ring and hence consider the affine cone
  dim(J) == 3 || error("not a surface")
  JM = jacobi_matrix(gens(J))
  minvec = minors(JM,ncharts - 3)
  I_sl = ideal(R,minvec) + J
  dim(I_sl) > 0 || return true                 ## smooth

  decomp = minimal_primes(I_sl)

  not_over_k = Vector{Tuple{MPolyIdeal,Int}}()
  dehom_dict = IdDict{Int,MatrixElem}()
  Jdehom_dict = IdDict{Int,MPolyIdeal}()

  for I in decomp
    d = dim(I)
    d < 2 || return false                      ## non-isolated is not duVal
    
    ## It is definitely a point now, i.e. I becomes maximal after dehomog.
    ## It shows up in at least one dehomogenization chart
    i = 0
    I_dehomog = dehomogenization(I,1)
    while i < ncharts
      i=i+1
      I_dehomog = dehomogenization(I,i)
      is_one(I_dehomog) || break
    end

    ## cache dehomogenized matrix
    if !haskey(dehom_dict,i)
      dehom_dict[i] = map(x -> dehomogenization(x,i), JM[:,:])
      Jdehom_dict[i] = dehomogenization(J,i)
    end

    ## make sure to do the geometrically irreducible points first
    if vdim(quo(base_ring(I_dehomog),I_dehomog)[1]) > 1
      push!(not_over_k, (I_dehomog,i))
      continue
    end

    _check_duval_at_point(I_dehomog,dehom_dict[i],Jdehom_dict[i]) || return false

  end

  ## handle geometrically reducible points afterwards
  for (I,i) in not_over_k
    li = absolute_primary_decomposition(I)
    for (_, _, I2, _) in li
      r_changed = base_ring(I2)
      kk = coefficient_ring(r_changed)
## Question: The following line is ugly, what is the nice way to move a matrix through change of coeff. field
      JM_changed = matrix(r_changed, ncharts, ngens(J), [change_coefficient_ring(kk,a) for a=dehom_dict[i][:,:]])
      J_changed = ideal(r_changed, [change_coefficient_ring(kk,a) for a=gens(Jdehom_dict[i])])
      _check_duval_at_point(I2,JM_changed,J_changed)
    end
  end

  return true
end

## CAUTION!!!
## The following needs to be cleaned up properly after vdim for modules is available in Oscar
##
function _check_duval_at_point(I_r::MPolyIdeal, JM_dehomog::MatrixElem, J_dehomog::MPolyIdeal)  
  ## go to chosen affine chart
  r = base_ring(I_r)
  I_max = ideal(r,gens(r))

  ## localize at point
  a = Oscar.rational_point_coordinates(I_r)
  U = complement_of_point_ideal(r,a)
  rl, loc_map = localization(r,U)
## Question: This following line is ugly, what is the nice way to move a matrix through extension of scalars
  JM_loc = map(x ->loc_map(x), JM_dehomog[:,:])

  if !all(iszero, a)
    Jm = Oscar.SubModuleOfFreeModule(JM_loc)
    F = ambient_free_module(Jm)
    Jm_shifted = Oscar.shifted_module(Jm)
    F_shifted = ambient_free_module(Jm_shifted)
    Jm_shifted = Jm_shifted + Oscar.shifted_ideal(loc_map(J_dehomog)) * F_shifted
  else
    Jm_shifted = Oscar.SubModuleOfFreeModule(JM_dehomog)
    F_shifted = ambient_free_module(Jm_shifted)
    Jm_shifted = Jm_shifted + J_dehomog * F_shifted
  end

  o = negdegrevlex(r)*lex(F_shifted)
  F1 = leading_module(Jm_shifted,o)
  F1_red = [reduce(gen(F_shifted,i),F1) for i in 1:ngens(F_shifted)]
  F1_red = filter(x -> !is_zero(x), F1_red)

  length(F1_red) == 1 || return false  ## not a hypersurface singularity
  
  F2 = [reduce(a*gen(F_shifted,i),F1) for i in 1:ngens(F_shifted) for a in gens(I_max)]
  F2_red = filter(x -> !is_zero(x),F2)
  corank = length(F2_red)
  # corank >=3 implies at least T_3,3,3, not duVal
  corank < 3 || return false   
  # corank ==1 implies A series        
  corank > 1 || return true    
   
  F3 = [reduce(a*gen(F_shifted,i),F1) for i in 1:ngens(F_shifted) for a in gens(I_max^2)]
  F3_red = filter(x -> !is_zero(x),F3)
  cubiccount = length(F3_red)
  # corank == 2 and cubiccount > 2 implies  at least X_9, not duVal 
  cubiccount < 3 || return false       
  # corank == 2 and cubiccount < 2 implies D series
  cubiccount > 1 || return true        

  # we are in E/J series
#  tau = vdim(SubquoModule(F,F1)[1])
  F4 = [reduce(a*gen(F_shifted,i),F1) for i in 1:ngens(F_shifted) for a in gens(I_max^3)]
  F4_red = filter(x -> !is_zero(x),F4)
  F5 = [reduce(a*gen(F_shifted,i),F1) for i in 1:ngens(F_shifted) for a in gens(I_max^4)]
  F5_red = filter(x -> !is_zero(x),F5)
  F6 = [reduce(a*gen(F_shifted,i),F1) for i in 1:ngens(F_shifted) for a in gens(I_max^5)]
  F5_red = filter(x -> !is_zero(x),F6)

  # tau > 8 implies at least J_10, not du Val   
  # worksaround for missing vdim: 
  length(F4_red) + length(F5_red) + length(F6_red) < (9-5) || return false             
  # E_k, k \in {6,7,8}
  return true                          
end

