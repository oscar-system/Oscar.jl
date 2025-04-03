###################################################################
## Tools for maximal contact objects                              
## - getters for Scheme, Covering, current ambient space etc
## - constructors  and updating 
###################################################################

###################################################################
## Getters
###################################################################
original_scheme(MCO::Oscar.MaxContactObject) = MCO.W_orig
covering(MCO::Oscar.MaxContactObject) = MCO.C
maximal_contact_data(MCO::Oscar.MaxContactObject) = MCO.max_contact_data

ambient_generators(MCC::Oscar.MaxContactChart) = MCC.ambient_gens
ambient_orders(MCC::Oscar.MaxContactChart) = MCC.ambient_orders
dependent_variables(MCC::Oscar.MaxContactChart) = MCC.dependent_vars
max_contact_minor_data(MCC::Oscar.MaxContactChart) = MCC.minor_data
prepared_jacobi_matrices(MCC::Oscar.MaxContactChart) = MCC.ambient_jacobi

function current_hypersurface_sequences(MC::Oscar.MaxContactObject)
  max_contact_data = maximal_contact_data(MC)
  res_dict = IdDict{AbsAffineScheme,Vector{MPolyRingElem}}()
  for U in keys(max_contact_data)
    res_dict[U] = ambient_generators(max_contact_data[U])
  end
  return res_dict
end

function current_order_sequences(MC::Oscar.MaxContactObject)
  max_contact_data = maximal_contact_data(MC)
  res_dict = IdDict{AbsAffineScheme, Vector{Int}}()
  for U in keys(max_contact_data)
    res_dict[U] = ambient_orders(max_contact_data[U])
  end
  return res_dict
end

###################################################################
## Constructors and Updating
###################################################################
function _initialize_max_contact_object(inc::Oscar.CoveredClosedEmbedding)

  W = codomain(inc)

  ## the empty hulls for the Oscar.MaxContactObject --  to be filled in
  max_contact_data = IdDict{AbsAffineScheme,Oscar.MaxContactChart}()
  patch_list = AbsAffineScheme[]

  ## run through the charts and set up a complete intersection covering
  ## for the locally complete intersection
  for U in affine_charts(W)
    IU = modulus(OO(U))

    if ngens(IU) == 0 
    ## trivial case -- nothing to be done
      max_contact_data[U] = Oscar.MaxContactChart(U,[],[],[],[],[])
      continue
    end

    ## find minors of jacobian matrix involved in expressing 1
    RU = base_ring(OO(U))
    IUgens = lifted_numerator.(gens(IU))
    JM = jacobi_matrix(IUgens)
    min_list = [a for a in minors_with_position(JM, nvars(RU) - dim(OO(U))) if !is_zero(a[1])]
    min_id = ideal(RU,[a[1] for a in min_list])
    min_id_full = min_id+ideal(RU,IUgens)+modulus(OO(U))
    radical_membership(one(RU),min_id_full)  || error("image of embedding not smooth")
    coord_vec = coordinates(one(RU),min_id_full)
    nonzero_indices =  [k for k in 1:ngens(min_id) if !is_zero(coord_vec[k])]
    
    ## for each refined chart fill in the data
    for k in nonzero_indices
      current_patch = PrincipalOpenSubset(U,OO(U).(min_list[k][1]))
      push!(patch_list, current_patch)
      amb_rows = min_list[k][2]
      amb_cols = min_list[k][3]
      selected_gens = [IUgens[i] for i in amb_cols]             ## selected generators via matrix columns
      if length(amb_cols) < ncols(JM)
        JM_essential = JM[:, amb_cols]
      else
        JM_essential = JM
      end
      submat_for_minor = JM[min_list[k][2], min_list[k][3]]
      Ainv, _ = pseudo_inv(submat_for_minor)
      JM_essential = JM_essential * Ainv          
      max_contact_data[current_patch] = Oscar.MaxContactChart(
                            U, selected_gens, 
                            [0 for i in 1:(nvars(RU) - dim(OO(U)))],
                            min_list[k][2],
                            [(min_list[k][1],length(min_list[k][2]))],
                            [JM_essential])
    end
  end
  
  new_Cov = Covering(patch_list)
  Oscar.inherit_gluings!(new_Cov,default_covering(W))
  return Oscar.MaxContactObject(W,new_Cov,max_contact_data)
end

function _max_contact_step(MC::Oscar.MaxContactObject,I::IdealSheaf,b::Int)
  Cov = covering(MC)
  max_contact_data = maximal_contact_data(MC)

  new_max_contact = IdDict{AbsAffineScheme,Oscar.MaxContactChart}()
  patch_list = AbsAffineScheme[]

  for U in Cov
    if is_zero(I(U)) || is_one(I(U))
      ## nothing to do on this chart 
      new_max_contact[U] = max_contact_data[U]
      continue
    end

    RU = base_ring(OO(U))
    IUgens = lifted_numerator.(gens(I(U)))
    JM = jacobian_matrix(max_contact_data[U], I(U))
    i = 1

    ## find max number of independent max contact hypersurfaces from I(U)
    save_list = []
    save_id = ideal(OO(U),[])
    ambient_contribution = ideal(RU,IUgens) + ideal(RU,lifted_numerator.(gens(modulus(OO(U)))))
    while i < min(nrows(JM),ncols(JM))
      ## find largest minors providing unit ideal
      testing_list = [a for a in minors_with_position(JM,i) if !is_zero(a[1])]
      testing_ideal =  ideal(RU, [a[1] for a in testing_list])
      radical_membership(_agnostic_complement_equation(U),
                             testing_ideal+ambient_contribution) || break 
      save_list = testing_list
      save_id = testing_ideal
      i = i+1
    end

    coord_vec = coordinates(_agnostic_complement_equation(U),
                             save_id+ambient_contribution)          # combination of 1
    nonzero_indices = [k for k in 1:ngens(save_id) if !is_zero(coord_vec[k])]  
                                                                    # indices of minors involved

    ## potential improvement: sort nonzero_indices by number of terms in minors
    cover_test = MPolyRingElem[]

    # use D(save_list[k][1]) for refinement -- k running through all non-zero indices
    for k in nonzero_indices
      push!(cover_test,save_list[k][1])

      ## refine U and update data
      if length(nonzero_indices) != 1
        current_patch = PrincipalOpenSubset(U,OO(U)(save_list[k][1]))
      else
        ## no refinement necessary, one minor is an invertible field element
        current_patch = U
      end

      push!(patch_list, current_patch)

      ## append new hypersurface generators and their orders
      MCU = max_contact_data[U]
      amb_gens = copy(ambient_generators(MCU))
      append!(amb_gens, [IUgens[a] for a in save_list[k][2]])
      amb_ord = copy(ambient_orders(MCU))
      append!(amb_ord, [b for a in save_list[k][2]])

      ## append row indices, new minor 
      dep_vars = copy(dependent_variables(MCU))
      append!(dep_vars,save_list[k][2])
      minor_data = copy(max_contact_minor_data(MCU))    ## this should save MCU.minor_data from getting altered?
      push!(minor_data,(save_list[k][1],length(save_list[k][2])))

      ## prepare and append relevant part of JM * pseudoinverse
      if length(save_list[k][3]) < length(IUgens)
        JM_essential = JM[:, save_list[k][3]]
      else
        JM_essential = JM
      end
      submat_for_minor = JM[save_list[k][2], save_list[k][3]]
      Ainv, _ = pseudo_inv(submat_for_minor)
      JM_essential = JM_essential * Ainv         # appropriate submatrix transformed to h*unit_matrix
      amb_jacobi = prepared_jacobi_matrices(MCU)
      push!(amb_jacobi, JM_essential)

      MCC_new = Oscar.MaxContactChart(U,amb_gens, amb_ord, dep_vars, minor_data, amb_jacobi)
      new_max_contact[current_patch] = MCC_new

      ## check whether we can stop prematurely, because everything is covered
      !is_empty(subscheme(U,cover_test)) || break
    end
  end
  
  
  ## create new MaxContactObject
  newCov = Covering(patch_list)
  return Oscar.MaxContactObject(original_scheme(MC),
                            newCov,
                            new_max_contact) 
end

_agnostic_complement_equation(U::PrincipalOpenSubset) = lift(complement_equation(U))
_agnostic_complement_equation(X::AbsAffineScheme{<:Field, <:MPolyRing}) = one(OO(X))
_agnostic_complement_equation(X::AbsAffineScheme) =  one(base_ring(OO(X)))
########################################################
## jacobian matrix in chart in maximal contact setting
########################################################
function jacobian_matrix(MCU::Oscar.MaxContactChart, I::Ideal)

  JI = jacobi_matrix(lifted_numerator.(gens(I)))
  var_indices = dependent_variables(MCU)
  minors_and_blocks = max_contact_minor_data(MCU)
  reduction_matrices = prepared_jacobi_matrices(MCU)
  
  length(minors_and_blocks) == length(reduction_matrices) || error("inconsistent ambient parameter data")
  cur_pos = 1
  old_mat = JI
  res_mat = old_mat
  for a in 1:length(minors_and_blocks)
  ## adjust to ambient space from maximal contact data
  ## minors and blocks[a] holds the minor and its size from the a-th step
  ## we need ot adjust for each step, one after the other
    res_mat = old_mat * minors_and_blocks[a][1]
    for i in 1:minors_and_blocks[a][2]
      for j in 1:ncols(res_mat)
    ## cur_pos holds the current total position, i the position in this block
    ## var_indices[cur_pos] holds the index of the current (dependent) variable to adjust with 
      res_mat[:,j] = res_mat[:,j] - 
         [reduction_matrices[a][k,i] * old_mat[var_indices[cur_pos],j] 
                     for k in 1:nrows(reduction_matrices[a])]
      end
      cur_pos = cur_pos + 1
    end
    old_mat = res_mat
  end

  return(res_mat)
end

function _delta_ideal_sheaf(MC::Oscar.MaxContactObject,I_sheaf::Oscar.AbsIdealSheaf)
  Cov = covering(MC)
  I_simple = small_generating_set(I_sheaf)

  ## run through all charts, compute the jacobian matrices w.r.t. the respective system of parameters
  ## and flatten the matrix to append it to gens(I)
  Delta_dict = IdDict{AbsAffineScheme,Ideal}()
  for U in Cov
    if is_one(I_simple(U))
      Delta_dict[U] = I_simple(U)
      continue
    end

    MCU = maximal_contact_data[U]
    result_mat = jacobian_matrix(MCU, I_simple(U))
    Ivec = copy(gens(I))
    append!(Ivec, ambient_generators(MCU))
    append!(Ivec,[a for a in collect(result_mat)])
    Delta_dict[U] = ideal(OO(U), OO(U).(Ivec))
  end

  return small_generating_set(IdealSheaf(scheme(I_sheaf),Delta_dict))
end

function _delta_list(MC::Oscar.MaxContactObject,I_sheaf::Oscar.AbsIdealSheaf,b::Int=0)
  Delta_list = AbsIdealSheaf[]
  j = 0
  I = I_sheaf
  while (( !is_one(I) && b == 0 ) || ( j < b ))
    push!(Delta_list, I)
    I = _delta_ideal_sheaf(MC,I)
    j = j + 1
  end

  return Delta_list
end
############################################################################################
## Maximal contact loci -- in different context
############################################################################################
function _nu_star_not_one_max(I::Oscar.AbsIdealSheaf)
  W = scheme(I)
  Cov = covering(W)
  inc_W = Oscar.CoveredClosedEmbedding(W,I)

  MC = _initialize_max_contact_object(inc_W)   ## only taking into account ambient W
  DI = _delta_ideal_sheaf(MC, I)
  if !is_one(DI)
    DL = _delta_list(MC.DI)
    return DL[:end],0,length(DL)+1
  end

  ## we can be sure that the maximal order of I is one
  MC2 = _max_contact_step(MC,I)               ## first maximal contact step

  shortest = dim(base_ring(OO(W)))            ## how many 1 in nu^*
  ret_dict = IdDict{AbsAffineScheme,Ideal}()

  for U in covering(MC2)
    MC2U = maximal_contact_data(MC2)[U]
    temp_len = length([a for a in ambient_orders(MC2U) if a == 1])     ## not counting leading zeros from W
    shortest >= temp_len || continue                                   ## not relevant for maximal value
    if shortest > temp_len
      shortest  = temp_len                                             ## new shortest one
      ret_dict = IdDict{AbsAffineScheme,Ideal}()
    end
    ret_dict[U] = ideal(OO(U),OO(U).(ambient_gens(MC2U)))
  end

  ideal_dict = IdDict{AbsAffineSchme,Ideal}()
  for U in covering(MC2)
    if U in keys(ret_dict)
      ideal_dict[U] = I(U)
    else
      I_one = ideal(OO(U), one(OO(U)))
      ret_dict[U] = I_one
      ideal_dict[U] = I_one
    end
  end

  J = IdealSheaf(W,red_dict;check=false)
  I_new = IdealSheaf(W,ideal_dict;check=false)
  DL = _delta_list(MC2,I_new)
  return J,DL[:end],shortest,length(DL),MC2
end
