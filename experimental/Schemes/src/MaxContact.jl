###################################################################
## Tools for maximal contact objects                              
## - getters for Scheme, Covering, current ambient space etc
## - constructors  and updating 
###################################################################

###################################################################
## Getters
###################################################################
original_scheme(MC::Oscar.MaxContactObject) = MC.W_orig
covering(MC::Oscar.MaxContactObject) = MC.C
maximal_contact_data(MC::Oscar.MaxContactObject) = MC.max_contact_data
ambient_parameter_data(MC::Oscar.MaxContactObject) = MC.ambient_param_data

function current_hypersurface_sequences(MC::Oscar.MaxContactObject)
  max_contact_data = maximal_contact_data(MC)
  res_dict = IdDict{AbsAffineScheme,Vector{MPolyRingElem}}()
  for U in keys(max_contact_data)
    res_dict[U] = max_contact_data[U][1]
  end
  return res_dict
end

function current_order_sequences(MC::Oscar.MaxContactObject)
  max_contact_data = maximal_contact_data(MC)
  res_dict = IdDict{AbsAffineScheme, Vector{Int}}()
  for U in keys(max_contact_data)
    res_dict[U] = max_contact_data[U][3]
  end
  return res_dict
end

function restricted_ideal_data(MC::Oscar.MaxContactObject)
  max_contact_data = maximal_contact_data(MC)
  res_dict = IdDict{AbsAffineScheme,Vector{MPolyRingElem}}()
  for U in keys(max_contact_data)
    res_dict[U] = max_contact_data[U][2]
  end
  return res_dict
end

###################################################################
## Constructors and Updating
###################################################################
function _initialize_max_contact_object(inc::Oscar.CoveredClosedEmbedding)

  W = codomain(inc)

  ## the empty hulls for the Oscar.MaxContactObject --  to be filled in
  max_contact_data = IdDict{AbsAffineScheme,
                            Tuple{Vector{MPolyRingElem},Vector{Int}}}()
  ambient_param_data = IdDict{AbsAffineScheme,
                            Tuple{Vector{Int},Vector{Tuple{RingElem,Int}},Vector{MatrixElem}}}()
  patch_list = AbsAffineScheme[]

  ## run through the charts and set up a complete intersection covering
  ## for the locally complete intersection
  for U in affine_charts(W)
    IU = modulus(OO(U))

    if ngens(IU) == 0 
    ## trivial case -- nothing to be done
      max_contact_data[U] = ([],[])
      ambient_param_data[U] = ([],[],[])
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
      max_contact_data[current_patch] = (selected_gens, [0 for i in 1:(nvars(RU) - dim(OO(U)))])
                                                                ## entry 2 marks that all are from orig W
      if length(amb_cols) < ncols(JM)
        JM_essential = JM[:, amb_cols]
      else
        JM_essential = JM
      end
      submat_for_minor = JM[min_list[k][2], min_list[k][3]]
      Ainv, _ = pseudo_inv(submat_for_minor)
      JM_essential = JM_essential * Ainv          
      ambient_param_data[current_patch] = (min_list[k][2], [(min_list[k][1],length(min_list[k][2]))], 
                                           [JM_essential])
    end
  end
  
  new_Cov = Covering(patch_list)
  Oscar.inherit_gluings!(new_Cov,default_covering(W))
  return Oscar.MaxContactObject(W,new_Cov,max_contact_data, ambient_param_data)
end

function _max_contact_step(MC::Oscar.MaxContactObject,I::IdealSheaf,b::Int)
  Cov = covering(MC)
  max_contact_data = maximal_contact_data(MC)
  ambient_param_data = ambient_parameter_data(MC)

  new_max_contact = IdDict{AbsAffineScheme,Tuple{Vector{MPolyRingElem},Vector{Int}}}()
  new_ambient_param_data = IdDict{AbsAffineScheme, 
                                  Tuple{Vector{Int},Vector{Tuple{RingElem,Int}},
                                  Vector{<:MatrixElem}}}()
  patch_list = AbsAffineScheme[]

  for U in Cov

    if is_zero(I(U)) || is_one(I(U))
      ## nothing to do on this chart 
      new_max_contact[U] = max_contact_data[U]
      new_ambient_param_data[U] = ambient_param_data[U]
      continue
    end

    RU = base_ring(OO(U))
    IUgens = lifted_numerator.(gens(I(U)))
    JM = jacobian_matrix(MC, I, U)
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
      temp_data = (copy(max_contact_data[U][1]),copy(max_contact_data[U][2]))
      append!(temp_data[1], [IUgens[a] for a in save_list[k][2]])
      append!(temp_data[2], [b for a in save_list[k][2]])
      new_max_contact[current_patch] = temp_data

      ## append row indices, new minor 
      temp_data = (copy(ambient_param_data[U][1]),copy(ambient_param_data[U][2]),copy(ambient_param_data[U][3]))
      append!(temp_data[1], save_list[k][2]) 
      push!(temp_data[2], (save_list[k][1],length(save_list[k][2])))

      ## prepare and append relevant part of JM * pseudoinverse
      if length(save_list[k][3]) < length(IUgens)
        JM_essential = JM[:, save_list[k][3]]
      else
        JM_essential = JM
      end
      submat_for_minor = JM[save_list[k][2], save_list[k][3]]
      Ainv, _ = pseudo_inv(submat_for_minor)
      JM_essential = JM_essential * Ainv         # appropriate submatrix transformed to h*unit_matrix
      push!(temp_data[3], JM_essential)

      new_ambient_param_data[current_patch] = temp_data

      ## check whether we can stop prematurely, because everything is covered
      !is_empty(subscheme(U,cover_test)) || break
    end
  end
  
  
  ## create new MaxContactObject
  newCov = Covering(patch_list)
  return Oscar.MaxContactObject(original_scheme(MC),
                            newCov,
                            new_max_contact,
                            new_ambient_param_data) 
end

_agnostic_complement_equation(U::PrincipalOpenSubset) = lift(complement_equation(U))
_agnostic_complement_equation(X::AbsAffineScheme{<:Field, <:MPolyRing}) = one(OO(X))
_agnostic_complement_equation(X::AbsAffineScheme) =  one(base_ring(OO(X)))
########################################################
## jacobian matrix in chart in maximal contact setting
########################################################
function jacobian_matrix(MC::Oscar.MaxContactObject, I::IdealSheaf, U::AbsAffineScheme)
  U in covering(MC) || error("U is not a chart of the given maximal contact covering")
  max_contact_data = maximal_contact_data(MC)[U]
  amb_data = ambient_parameter_data(MC)[U]

  JI = jacobi_matrix(lifted_numerator.(gens(I(U))))
  var_indices = amb_data[1]
  minors_and_blocks = amb_data[2]
  reduction_matrices = amb_data[3]
  
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