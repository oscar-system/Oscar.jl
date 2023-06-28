function has_du_val_singularities(X::AbsProjectiveScheme{<:Field,<:Any})
  return has_du_val_singularities(covered_scheme(X))
end

function has_du_val_singularities(X::AbsSpec{<:Field,<:Any})
  R = OO(X)
  I = modulus(R)
  
  J = image_ideal(singular_locus(X)[2])
  dim(J) == 0 || return false                            ## non-isolated
  return is_du_val_singularity(X,J)
end

function has_du_val_singularities(X::AbsCoveredScheme{<:Field})
  C = (has_attribute(X0, :simplified_covering) ? simplified_covering(X0) : default_covering(X0))

  I = ideal_sheaf_of_singular_locus(X)
  decomp = minimal_associated_points(I)

  ## we do the double loop here to avoid unnecessary checks
  for J in decomp
    for U in C
      !isone(J(U)) || continue
      is_du_val_singularity(X(U),J(U)) || return false    ## testing the point in one chart suffices
      break                         
    end
  end

  return true
end

function is_du_val_singularity(X::AbsSpec{<:Field,<:Any},I::Union{MPolyIdeal, MPolyLocalizedIdeal})
  OOX = OO(X)
  dim(X) == 2 || error("Scheme not of dimension 2")
  J = modulus(OOX)
  !isone(J) || error("Scheme is empty")
  !iszero(J) || return true                  ## X smooth

  R = base_ring(I)
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")
  base_ring(OOX) === R || error("base rings need to be identical")

  dim(I) == 0 || error("second argument is not a point")
  if get_attribute(I,:is_absolutely_prime, false)
    return _check_duval_at_point(J,I)[1]
  end

  decomp = absolute_primary_decomposition(I)
  for (_,_,I2,mult) in decomp
    set_attribute!(I,:is_absolutely_prime,true)
    ## pass to algebraic extension
    r_changed = base_ring(I2)
    kk = coefficient_ring(r_changed)
    J_changed = ideal(r_changed,  [change_coefficient_ring(kk,a, parent = r_changed) for a=gens(J)])
    is_du_val_singularity(Spec(quo(r_changed,J_changed)[1]),I2) || return false
  end
  
  return true
end

function decide_duval_at_point(X::AbsSpec{<:Field,<:Any},I::Union{MPolyIdeal, MPolyLocalizedIdeal})
  OOX = OO(X)
  dim(X) == 2 || error("Scheme not of dimension 2")
  J = modulus(OOX)
  !isone(J) || error("Scheme is empty")
  !iszero(J) || return true                  ## X smooth

  R = base_ring(I)
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")
  base_ring(OOX) === R || error("base rings need to be identical")

  dim(I) == 0 || error("second argument is not a point")
  if get_attribute(I,:is_absolutely_prime, false)
    li = _check_duval_at_point(J,I)
    return [(li[1], I, li[2])]
  end

  result_vector = []
  decomp = absolute_primary_decomposition(I)
  for (_,_,I2,mult) in decomp
    set_attribute!(I,:is_absolutely_prime,true)
    ## pass to algebraic extension
    r_changed = base_ring(I2)
    kk = coefficient_ring(r_changed)
    J_changed = ideal(r_changed,  [change_coefficient_ring(kk,a, parent = r_changed) for a=gens(J)])
    tempvec = decide_duval_at_point(Spec(quo(r_changed,J_changed)[1]),I2)
    for x in tempvec
      x[1] || return [x]
      push!(result_vector,x)
    end
  end
@show result_vector

  return result_vector
end

function _check_duval_at_point(IX::Ideal,Ipt::Ideal)
  R = base_ring(IX)
  R == base_ring(Ipt) || error("basering mismatch")
  kk = base_ring(R)
  characteristic(kk) == 0 || error("only available in characteristic zero")
 
  JM = jacobi_matrix(gens(IX))

  smooth = (:A,0)

  ## localize at point
  a = Oscar.rational_point_coordinates(Ipt)
  U = complement_of_point_ideal(R,a)
  RL, loc_map = localization(R,U)
  IX_loc = loc_map(IX)
  JM_loc =  map(x ->loc_map(x), JM[:,:])

  if !all(iszero(a))
    F_loc = free_module(RL,ngens(IX))
    Jm = sub(F_loc,JM_loc)
    Jm = Jm + loc_map(IX)*F_loc
    Jm_shifted = Oscar. shifted_module(Jm)
    F_shifted = ambient_free_module(Jm_shifted)
  else
    F = free_module(R,ngens(IX))
    Jm = sub(F,JM)[1]
    Jm_shifted = Jm + (IX * F)[1]
    F_shifted = F
  end

  o = negdegrevlex(R)*lex(F_shifted)
  F1 = leading_module(Jm_shifted,o)
  F1quo = quo(F_shifted, F1)[1]

  constant_mons = vector_space_dimension(F1quo,0) 
  constant_mons < 2 || return (false, typeof(smooth))                   ## not a hypersurface
  constant_mons > 0 || return (true, smooth)                            ## no singularity
  
  tau = vector_space_dimension(F1quo)

  corank = vector_space_dimension(F1quo,1)
  corank < 3 || return (false,typeof(smooth))                           ## at least T_3,3,3 not duVal
  corank > 1 || return (true,(:A,tau))                                  ## A_series

  # we now already know essentially a hypersurface of corank 2, count degrees of freedom cubicpart
  cubiccount = vector_space_dimension(F1quo,2)
  cubiccount < 3 || return  (false, typof(smooth))                      ## at least X_9
  cubiccount > 1 || return  (true, (:D,tau))                            ## D_series

  # it is definitely in the E/J series
  tau < 9 || return(false, typeof(smooth))                          ## at least J_10
  return (true, (:E,tau))                                               ## E_6, E_7, E_8
end

function vector_space_dimension(M::SubquoModule)
  
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  has_monomials_on_all_axes(LM) || error("not a finite dimensional vector space")
  
  d = 0
  sum_dim = 0
  tempdim = vector_space_dimension(M,0)

  while tempdim > 0
    sum_dim = sum_dim + tempdim
    d = d+1
    tempdim = vector_space_dimension(M,d)
  end
 
  return sum_dim
end

function vector_space_dimension(M::SubquoModule,d::Int64)
  R = base_ring(M)
  F = ambient_free_module(M)
  Mq,_ = sub(F,rels(M))

  ambient_representatives_generators(M) == gens(F) || error("not implemented for M/N with non-trivial M")

  o = default_ordering(M)
  LM = leading_module(Mq,o)

  return length([x*e for x in Oscar.all_monomials(R, d) for e in gens(F) if !(x*e in LM)])
end
  
function has_monomials_on_all_axes(M::SubquoModule)
  R = base_ring(M)

  length(rels(M)) == 0 || error("not implemented for quotients")
  
  ambient_rank = ngens(ambient_free_module(M))
  genlist = ambient_representatives_generators(M)
  explist = Tuple{Vector{Int64}, Int64, Int}[]
  for x in genlist
    tempexp = leading_exponent(x)
    tempdeg = sum(tempexp[1])
    push!(explist,(tempexp[1],tempexp[2],tempdeg))
  end
  all(findfirst(x -> (x[1][i] == x[3] && x[2]==j),explist) != nothing for i in 1:ngens(R) for j in 1:ambient_rank) || return false
  return true
end

