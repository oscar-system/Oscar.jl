## Warnung: show auf desingMor geht noch nicht!!!
export _desing_curve
export  find_refinement_with_local_system_of_params

#####################################################################################################
# Desingularization morphism: birational map between covered schemes with smooth domain
#####################################################################################################

# Fehlt: NormalizationMorphism fuer Schemata -- muessten wir haben, sobald wir Lipman machen wollen
#
#@attributes mutable struct LipmanStyleSequence{
#    DomainType<:AbsCoveredScheme,
#    CodomainType<:AbsCoveredScheme
#   } <: AbsDesingMor{
#                                 DomainType,
#                                 CodomainType,
#                    }
#  maps::Vector{:<union{BlowupMorphism,NormalizationMorphism}}        # count right to left:
#                                                 # original scheme is codomain of map 1
#  # boolean flags
#  resolves_sing::Bool                            # domain not smooth yet?
#  is_trivial::Bool                               # codomain already smooth?
#
#  # fields for caching, may be filled during computation
#  ex_div::Vector{<:EffectiveCartierDivisor}      # list of exc. divisors arising from individual steps
                                                  # in domain(maps[end])
#
#  # fields for caching to be filled a posteriori (on demand, only if partial_res==false)
#  composed_map::AbsCoveredSchemeMorphism        
#  exceptional_divisor::WeilDivisor          
#
#  function LipmanStyleSequence(maps::Vector{<:AbsCoveredSchemeMorphism})
#    n = length(maps)
#    for i in 1:n-1
#      @assert domain(maps[i]) === codomain(maps[i+1]) "not a sequence of morphisms"
#    end
#    return new{typeof(domain(maps[end])),typeof(codomain(first(maps)))}(maps)
#  end
#end


##################################################################################################
# getters
##################################################################################################
morphisms(phi::AbsDesingMor) = copy(phi.maps)
morphism(phi::AbsDesingMor,i::Int) = copy(phi.maps[i])
last_map(phi::AbsDesingMor) = phi.maps[end]
exceptional_divisor_list(phi::BlowUpSequence) = phi.ex_div  ## derzeit Liste von Eff. Cartier Div.
embeddings(phi::BlowUpSequence) = phi.embeddings

## do not use!!! (for forwarding and certain emergenies)
function underlying_morphism(phi::AbsDesingMor)
  if !isdefined(phi, :underlying_morphism)
    phi.underlying_morphism = CompositeCoveredSchemeMorphism(reverse(morphisms(phi)))
  end
  return phi.underlying_morphism
end

##################################################################################################
# setting values in DesingMors -- Watch out: only place with direct access to fields!!!
##################################################################################################
function add_map!(f::BlowUpSequence, phi::BlowupMorphism)
  push!(f.maps, phi)
  ex_div = [strict_transform(phi,E) for E in f.ex_div[1:end]]
  push!(ex_div, exceptional_divisor(phi))
  f.ex_div = ex_div
  return f
end

function initialize_blow_up_sequence(phi::BlowupMorphism)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_trivial = is_one(center(phi))
  f.resolves_sing = false                                # we have no information, wether we are done
                                                         # without further computation
  f.is_embedded = false
  return f
end

function add_map_embedded!(f::BlowUpSequence, phi::BlowupMorphism)
  push!(f.maps, phi)
  ex_div = typeof(f.ex_div[1])[strict_transform(phi, E) for E in f.ex_div[1:end]]
  push!(ex_div, exceptional_divisor(phi))
  f.ex_div = ex_div
  if f.transform_type == :strict
    X_strict, inc_strict,_ = strict_transform(phi, f.embeddings[end])
    push!(f.embeddings, inc_strict)
  elseif f.transform_type == :weak
    I_trans,b = weak_transform_with_multiplicity(phi, f.controlled_transform)
    push!(f.ex_mult,b)
    f.controlled_transform = I_trans
  else
    I_trans = controlled_transform(phi, f.controlled_transform, f.ex_mult[end])
    f.controlled_transform = I_trans
    push!(f.ex_mult, f.ex_mult[end])
  end
  return f
end

function initialize_embedded_blowup_sequence(phi::BlowupMorphism, inc::CoveredClosedEmbedding)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_embedded = true
  f.transform_type = :strict
  if !is_one(center(phi))
    f.is_trivial = false
    X_strict,inc_strict,_ = strict_transform(phi,inc)
    f.embeddings = [f, inc_strict]
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
  else
    f.is_trivial = true
    f.embeddings = [inc, inc]
    f.resolves_sing = false
  end
  return f
end

function initialize_embedded_blowup_sequence(phi::BlowupMorphism, I::AbsIdealSheaf, b::Int)
  f = BlowUpSequence([phi])
  f.ex_div = [exceptional_divisor(phi)]
  f.is_embedded = true
  if !is_one(center(phi))
    f.is_trivial = false
    if b == 0
      I_trans, b = weak_transform_with_multiplicity(phi,I)
      f.transform_type = :weak
    elseif b > 0
      I_trans = controlled_transform(phi, I, b)
      f.transform_type = :controlled
    end
    f.controlled_transform = I_trans                     # CAUTION: b is considered set once and for all
    f.ex_mult = [b]
    f.resolves_sing = false                              # we have no information, whether we are done
                                                         # without further computation
  else
    f.is_trivial = true
    f.controlled_transform = I
    f.transform_type = :weak
    f.ex_mult = [0]
    f.resolves_sing = false
  end
  return f
end


##################################################################################################
# desingularization workers
##################################################################################################
function embedded_desingularization(f::Oscar.CoveredClosedEmbedding; algorithm::Symbol=:BEV)
  I_sl = ideal_sheaf_of_singular_locus(domain(f))

  ## trivial case: domain(f) was already smooth
  if is_one(I_sl)
    id_W = identity_blow_up(codomain(f))
    phi = initialize_embedded_blowup_sequence(id_W,f)
    phi.resolves_sing = true
    return phi
  end

  ## I_sl non-empty, we need to do something
  dimX = dim(domain(f))
  if dimX == 1
    return _desing_emb_curve(f,I_sl)
#  elseif ((dimX == 2) && (algorithm == :CJS))
#    return _desing_CJS(f)
#  elseif (algorithm == :BEV)
#    return _desing_BEV(f)
  end
# here the keyword algorithm ensures that the desired method is called
  error("not implemented yet")
end

function embedded_desingularization(inc::ClosedEmbedding; algorithm::Symbol=:BEV)
  return embedded_desingularization(CoveredClosedEmbedding(inc); algorithm)
end

function CoveredClosedEmbedding(inc::ClosedEmbedding; domain=CoveredScheme(domain(inc)), codomain=CoveredScheme(codomain(inc)))
  mor_dict = IdDict{AbsAffineScheme, ClosedEmbedding}(domain[1][1] => inc)
  cov_mor = CoveringMorphism(default_covering(domain), default_covering(codomain), mor_dict; check=false)
  return CoveredClosedEmbedding(domain, codomain, cov_mor; check=false)
end

function desingularization(X::AbsCoveredScheme; algorithm::Symbol=:Lipman)
  I_sl = ideal_sheaf_of_singular_locus(X)
  
  ## trivial case: X is already smooth
  if is_one(I_sl)
    id_X = identity_blow_up(X)
    maps = [id_X] 
    return_value = BlowUpSequence(maps)
    return_value.resolves_sing = true
    return_value.is_trivial = true
    return return_value
  end

  ## I_sl non-empty, we need to do something 
# here the keyword algorithm ensures that the desired method is called
  dimX = dim(X)
  if dimX == 1
    return _desing_curve(X, I_sl)
  end
#  if ((dimX == 2) && (algorithm==:Lipman))
#    error("not implemented yet")
#    return_value = _desing_lipman(X, I_sl)
#    return return_value
#  end
#  if ((dimX == 2) && (algorithm==:Jung))
#    error("not implemented yet")
#    return_value = _desing_jung(X)
#   end       
  error("not implemented yet")    
end

function desingularization(X::AbsAffineScheme; algorithm::Symbol=:BEV)
  return desingularization(CoveredScheme(X); algorithm)
end

function _desing_curve(X::AbsCoveredScheme, I_sl::AbsIdealSheaf)
  ## note: I_sl not unit_ideal_sheaf, because this has been caught before in desingularization(X) 
  decomp = maximal_associated_points(I_sl)
  I = small_generating_set(pop!(decomp))
  current_blow_up = blow_up(I)
  phi = initialize_blow_up_sequence(current_blow_up)::BlowUpSequence
  decomp = [strict_transform(current_blow_up,J) for J in decomp]
  
  I_sl_temp = I_sl
  while !is_one(I_sl_temp)
    while length(decomp) > 0
      I = small_generating_set(pop!(decomp))
      phi = _do_blow_up!(phi,I)
      if length(decomp)>0 
        decomp = [strict_transform(last_map(phi),J) for J in decomp]
      end
    end
    I_sl_temp = ideal_sheaf_of_singular_locus(domain(last_map(phi)))
    decomp = maximal_associated_points(I_sl_temp)
  end

  phi.resolves_sing = true
  return phi
end

function _desing_emb_curve(f::CoveredClosedEmbedding, I_sl::AbsIdealSheaf)
  ## note: I_sl not unit_ideal_sheaf, because this has been caught before in embedded_desingularization(f)
  decomp = maximal_associated_points(pushforward(f)(I_sl))
  I = small_generating_set(pop!(decomp))
  current_blow_up = blow_up(I)
  phi = initialize_embedded_blowup_sequence(current_blow_up,f)::BlowUpSequence
  decomp = [strict_transform(current_blow_up,J) for J in decomp]
  I_sl_temp = I_sl
  while !is_one(I_sl_temp)
    while length(decomp) > 0
      I = small_generating_set(pop!(decomp))
      phi = _do_blow_up_embedded!(phi,I)
      if length(decomp)>0
        decomp = [strict_transform(last_map(phi),J) for J in decomp]
      end
    end
    last_emb = embeddings(phi)[end]
    I_sl_temp = pushforward(last_emb, ideal_sheaf_of_singular_locus(domain(last_emb)))
    decomp = maximal_associated_points(I_sl_temp)
  end

  phi = _ensure_ncr!(phi)
  phi.resolves_sing = true
  return phi
end

function _ensure_ncr!(f::AbsDesingMor)
  current_divs = exceptional_divisor_list(f)

# this first step can only come into play, if the centers were not determined algorithmically
# it is ensured by all standard desingularization algorithms
  I_bad = non_snc_locus(current_divs)
  while !is_one(I_bad)
    decomp = maximal_associated_points(I_bad)
    while length(decomp)>0
      I = small_generating_set(pop!(decomp))
      f =_do_blow_up_embedded!(f,I)
      if length(decomp)>0
        decomp = [strict_transform(last_map(f),J) for J in decomp]
      end
    end
    I_bad = non_snc_locus(exceptional_divisor_list(f))
  end

# now ensure the transversal intersections with the strict transform
  current_divs = copy(exceptional_divisor_list(f))
  I_X = image_ideal(f.embeddings[end])
  while !is_empty(current_divs)
    one_div = popfirst!(current_divs)       ## we need a FIFO here, not the usual LIFO
    I_temp = I_X + ideal_sheaf(one_div)
    last_emb = embeddings(f)[end]
    inc_temp = CoveredClosedEmbedding(scheme(I_temp), I_temp)
    next_locus = ideal_sheaf_of_singular_locus(domain(inc_temp))
    decomp = maximal_associated_points(pushforward(inc_temp, next_locus ))
    while !is_empty(decomp)
      I = small_generating_set(pop!(decomp))
      f =_do_blow_up_embedded!(f,I)
      I_X = image_ideal(f.embeddings[end])
      current_divs = [strict_transform(last_map(f),J) for J in current_divs]
      push!(current_divs, exceptional_divisor_list(f)[end])
    end
  end

# finally make sure not too many exceptional divisors meet the strict transform in the same point
  n_max = dim(I_X)
  current_divs = copy(exceptional_divisor_list(f))
  _,inter_div = divisor_intersections_with_X(current_divs,I_X)
  while !is_empty(inter_div)
    cent = small_generating_set(pop!(inter_div))
    f =_do_blow_up_embedded!(f,cent)
    if length(inter_div)>0
      inter_div = [strict_transform(last_map(f),J) for J in inter_div]
    end
  end
  return f
end


function _do_blow_up!(f::AbsDesingMor, cent::AbsIdealSheaf)
  old_sequence = morphisms(f)
  X = domain(old_sequence[end])
  X === scheme(cent) || error("center needs to be defined on same scheme")
  current_blow_up = blow_up(cent,var_name=string("v", length(old_sequence), "_"))
  add_map!(f, current_blow_up)
  return(f)
end

function _do_blow_up_embedded!(phi::AbsDesingMor,cent::AbsIdealSheaf)
  old_sequence = morphisms(phi)
  X = domain(old_sequence[end])
  X === scheme(cent) || error("center needs to be defined on same scheme")
  current_blow_up = blow_up(cent,var_name=string("v", length(old_sequence), "_"))
  add_map_embedded!(phi, current_blow_up)
  return(phi)
end


###################################################################################################
# Should go to IdealSheaf.jl, when PR is ready to merge
###################################################################################################

function unit_ideal_sheaf(X::AbsCoveredScheme)
  dd = IdDict{AbsAffineScheme, Ideal}(U=>ideal(OO(U), [one(OO(U))]) for U in affine_charts(X))
  return IdealSheaf(X, dd, check=false)
end

function zero_ideal_sheaf(X::AbsCoveredScheme)
  dd = IdDict{AbsAffineScheme, Ideal}(U=>ideal(OO(U), elem_type(OO(U))[]) for U in affine_charts(X))
  return IdealSheaf(X, dd, check=false)
end

function identity_blow_up(X::AbsCoveredScheme)
  f = BlowupMorphism(X, unit_ideal_sheaf(X))
  return f
end

########################################################################
# Refinements to find local systems of parameters
########################################################################

function find_refinement_with_local_system_of_params(W::AbsAffineScheme{<:Field, <:MPolyRing}; check::Bool=true)
  U = PrincipalOpenSubset(W, one(OO(W)))
  res_cov = Covering([U])
  R = ambient_coordinate_ring(W)
  minor_dict = IdDict{AbsAffineScheme, Tuple{Vector{Int}, Vector{Int}, elem_type(R)}}()
  minor_dict[U] = (Int[], Int[], one(R))
  return res_cov, minor_dict
end

function find_refinement_with_local_system_of_params(W::AbsAffineScheme; check::Bool=true)
  @check is_smooth(W) "scheme must be smooth"
  @check is_equidimensional(W) "scheme must be equidimensional"
  mod_gens = lifted_numerator.(gens(modulus(OO(W))))::Vector{<:MPolyRingElem}
  # We run into difficulties for the zero ideal as a modulus.
  # To get the matrix of minors of jac(I) we use `induced_map_on_exterior_power` below.
  # It is mathematical convention that ⋀⁰R⁰ = R¹. But that's unfortunately not 
  # coherent with the generic indexing we use in the implementation below. 
  # Should this assertion lead to problems, one can still replace throwing an error 
  # by inserting the shortcut from the method above and returning that. 
  @assert !isempty(mod_gens) "method not implemented for empty modulus; try to create the scheme without modulus instead"
  R = ambient_coordinate_ring(W)
  M = jacobian_matrix(R, mod_gens)

  n = nrows(M) # the number of variables in the ambient_ring
  r = ncols(M) # the number of generators

  Rn = FreeMod(R, n)
  Rr = FreeMod(R, r)
  phi = hom(Rn, Rr, M)
  codim = n - dim(W)
  phi_cod = induced_map_on_exterior_power(phi, codim)
  M_ext = matrix(phi_cod)
  n_cod = nrows(M_ext)
  r_cod = ncols(M_ext)

  all_entries = Vector{Int}[[i, j] for i in 1:n_cod for j in 1:r_cod]
  M_ext_vec = elem_type(R)[M[i, j] for i in 1:n_cod for j in 1:r_cod]
  min_id = ideal(OO(W), M_ext_vec)
  lambda_vec = coordinates(one(OO(W)), min_id)
  lambda = elem_type(OO(W))[lambda_vec[(i-1)*r_cod + j] for i in 1:n_cod, j in 1:r_cod]
  
  nonzero_indices_linear = [k for k in 1:length(lambda_vec) if !is_zero(lambda_vec[k])]
  non_zero_indices = [[i, j] for i in 1:n_cod, j in 1:r_cod if !is_zero(lambda_vec[(i-1)*r_cod + j])]

  ref_patches = AbsAffineScheme[]
  minor_dict = IdDict{AbsAffineScheme, Tuple{Vector{Int}, Vector{Int}, elem_type(R)}}()
  for (i, j) in non_zero_indices
    h_ij = M_ext[i, j]
    U_ij = PrincipalOpenSubset(W, OO(W)(h_ij))
    I = ordered_multi_index(i, codim, n)
    J = ordered_multi_index(j, codim, r)
    push!(ref_patches, U_ij)
    minor_dict[U_ij] = (indices(I), indices(J), M_ext[i, j])
  end
  res_cov = Covering(ref_patches)
  inherit_glueings!(res_cov, Covering(W))
  return res_cov, minor_dict
end

function find_refinement_with_local_system_of_params_rec(
    W::AbsAffineScheme, 
    mod_gens::Vector{PolyType} = lifted_numerator.(gens(modulus(OO(W)))),
    row_ind::Vector{Int} = Int[],
    col_ind::Vector{Int} = Int[],
    trans_mat::MatrixElem{RingElemType} = change_base_ring(OO(W), jacobi_matrix(mod_gens));
    check::Bool=true
  ) where {PolyType <: MPolyRingElem, RingElemType <: RingElem}

  # End of recursion
  n = dim(ambient_coordinate_ring(W))
  if length(row_ind) == n - dim(W)
    return [(W, row_ind, col_ind, prod(trans_mat[row_ind[k], col_ind[k]] for k in 1:dim(W); init=one(OO(W))))]
  end

  # generate the unit ideal of OO(W) with the entries of trans_mat
  n = nrows(trans_mat)
  r = ncols(trans_mat)
  all_entries_ind = [[i, j] for i in 1:n if !(i in row_ind) for j in 1:r if !(j in col_ind)]
  all_entries = elem_type(OO(W))[trans_mat[i, j] for (i, j) in all_entries_ind]
  entry_id = ideal(OO(W), all_entries)
  lambda = coordinates(one(OO(W)), entry_id)

  non_zero_entries = [k for k in 1:length(lambda) if !is_zero(lambda[k])]

  loc_results = Tuple{<:AbsAffineScheme, Vector{Int}, Vector{Int}, <:RingElem}[]
  for k in non_zero_entries
    i, j = all_entries_ind[k]
    h_ij = trans_mat[i, j]
    U_ij = PrincipalOpenSubset(W, OO(W)(h_ij))
    res_mat = change_base_ring(OO(U_ij), trans_mat) # TODO: Avoid checks here
    new_row_ind = vcat(row_ind, [i])
    new_col_ind = vcat(col_ind, [j])
    
    # Do Gaussian elimination on the matrix to kill off the other entries in this row
    u = res_mat[i, j]
    inv_u = inv(u)
    for l in 1:r
      l in new_col_ind && continue
      res_mat = add_column!(res_mat, -inv_u * res_mat[i, l], j, l)
    end
    loc_results = vcat(loc_results, 
          find_refinement_with_local_system_of_params_rec(
              U_ij, mod_gens, new_row_ind, new_col_ind, res_mat; check=check
             )
         )
  end
  return loc_results
end

##################################################################################################
#  locus of order at least b and of maximal order
##################################################################################################

function max_order_locus(I::AbsIdealSheaf)
  return _delta_list(I)[end]
end

function locus_of_order_geq_b(I::AbsIdealSheaf, b::Int)
  return _delta_list(I,b)[end]
end

function _delta_ideal_for_order(inc::CoveredClosedEmbedding)

  W = codomain(inc)
  I = image_ideal(inc)

  Delta_dict = IdDict{AbsAffineScheme,Ideal}()

  for U in affine_charts(W)
    XU, inc_U = sub(U,I(U))
    Cov,Chart_dict = find_refinement_with_local_system_of_params(U)
    Delta_dict[U] = (_delta_ideal_for_order(CoveredClosedEmbedding(inc_U), Cov, Chart_dict))(U)
  end

  return IdealSheaf(W,Delta_dict;check=false)
end

function _delta_list(inc::CoveredClosedEmbedding)
  I = image_ideal(inc)
  return _delta_list(I)
end

function _delta_list(I::AbsIdealSheaf, b::Int=0)
  W = scheme(I)
  is_smooth(W) || error("ambient scheme needs to be smooth")
  Delta_list = AbsIdealSheaf[]
  j = 0
  while (( !is_one(I) && b == 0 ) || ( j < b ))
    push!(Delta_list, I)
    inc = CoveredClosedEmbedding(scheme(I),I)
    I = _delta_ideal_for_order(inc)
    j = j + 1
  end

  return Delta_list
end

function _delta_ideal_for_order(inc::CoveredClosedEmbedding, Cov::Covering, 
       ambient_param_data::IdDict{<:AbsAffineScheme,
                                 <: Tuple{Vector{Int64},Vector{Int64},<:RingElem}};
       check::Bool=true)

  W = codomain(inc)                                
  @check is_smooth(W) "codomain of embedding needs to be smooth"
#  @check is_equidimensional(W) "codomain of embedding needs to be equidimensional"
  I_X = small_generating_set(image_ideal(inc))         # ideal sheaf describing X on W

  Delta_dict = IdDict{AbsAffineScheme,Ideal}()
  for U in Cov
    I = I_X(U)
    if is_one(I)
      Delta_dict[U] = I
      continue
    end

    amb_row,amb_col,h = ambient_param_data[U]
    mod_gens = lifted_numerator.(gens(modulus(OO(U))))
    R = ambient_coordinate_ring(U)
    JM = jacobian_matrix(R, mod_gens)
    if length(amb_col) < length(mod_gens)
      JM_essential = JM[:, amb_col]
    else
      JM_essential = JM
    end
    submat_for_minor = JM[amb_row, amb_col]
    Ainv, h2 = pseudo_inv(submat_for_minor)
    h == h2 || error("inconsistent input data")
    JM_essential = JM_essential * Ainv
    I_gens = lifted_numerator.(gens(I))
    JI = jacobian_matrix(I_gens)
    result_mat = h*JI
    for i in 1:length(amb_col)
      for j in 1:ncols(result_mat)
        result_mat[:,j] = result_mat[:,j] - [JM_essential[k,i] * JI[amb_row[i],j] for k in 1:nrows(JM_essential)]
      end
    end
    Ivec = copy(gens(I))
    append!(Ivec,[a for a in OO(U).(collect(result_mat))])
    Delta_dict[U] = ideal(OO(U),Ivec)
  end

  return small_generating_set(IdealSheaf(W,Delta_dict))
end
 
########################################################################
# test for snc                                                         #
########################################################################

function divisor_intersections_with_X(current_div, I_X)
  scheme(I_X) == scheme(current_div[1]) || error("underlying schemes do not match")
  n_max = dim(I_X)

  inter_div_dict = Dict{Vector{Int},Tuple{AbsIdealSheaf,Int}}()
  old_keys = Vector{Int}[]
  empty_keys = Vector{Int}[]
  essential_inter = AbsIdealSheaf[]

# initialization: each divisor + I_X
  for k in 1:length(current_div)
    Idiv = ideal_sheaf(current_div[k]) + I_X
    if !is_one(Idiv)
      inter_div_dict[[k]] = (Idiv,0)
      push!(old_keys, [k])
    end
  end
  new_keys = copy(empty_keys)

# add intersections
  while !is_empty(old_keys)
    old_keyvec = popfirst!(old_keys)       # this is the intersection to which we add a new divisor
    for i in (old_keyvec[end]+1):length(current_div)
      haskey(inter_div_dict, [i]) || continue     # divisor i meets X at all
      copykey = copy(old_keyvec)
      push!(copykey,i)                     # here we add it
      Idiv = inter_div_dict[[i]][1] + inter_div_dict[old_keyvec][1]
      !is_one(Idiv) || continue            # it intersects

      subsetlist = subsets(old_keyvec,length(old_keyvec)-1)
      subsetlist = [push!(a,i) for a in subsetlist]
      if (sum([inter_div_dict[a][2] for a in subsetlist])> 0)
                                           # offending intersection, known before
        inter_div_dict[copykey] = (Idiv,2)
      elseif dim(Idiv) == n_max - length(copykey)
                                           # non-offending intersection
        inter_div_dict[copykey] = (Idiv,0)
      else
                                           # offending intersection, new
        inter_div_dict[copykey] = (Idiv,1)
        push!(essential_inter, Idiv)
      end
      push!(new_keys,copykey)
    end

    # go to next higher number of intersecting divisors
    if is_empty(old_keys)
      old_keys = copy(new_keys)
      new_keys = copy(empty_keys)
    end
  end

  return inter_div_dict, essential_inter
end

function non_snc_locus(divs::Vector{<:EffectiveCartierDivisor})
  is_empty(divs) && error("list of divisors must not be empty")
  X = scheme(first(divs))
  @assert all(d->scheme(d) === X, divs)
  @assert is_smooth(X)
  r = length(divs)
  triv_cov = trivializing_covering.(divs)

  com_ref, incs = common_refinement(triv_cov, default_covering(X))

  ideal_dict = IdDict{AbsAffineScheme, Ideal}() # ideal sheaf of the non_snc_locus
  for U in patches(com_ref)
    loc_eqns = elem_type(OO(U))[]
    for k in 1:length(incs)
      I = ideal_sheaf(divs[k])
      inc = incs[k][U]
      V = codomain(inc)
      h = gen(I(V),1)
      hh = pullback(inc)(h)
      is_unit(hh) && continue # Not every divisor needs to be visible here
      push!(loc_eqns, hh)
    end
    if isempty(loc_eqns) || is_regular_sequence(loc_eqns)
      ideal_dict[U] = ideal(OO(U), one(OO(U))) # Nothing visible here
      continue
    end

    # Determine the non-snc-locus
    K = koszul_complex(KoszulComplex, loc_eqns)
    k = findfirst(k->!is_zero(homology(K, k)[1]), 1:length(loc_eqns))
    ideal_dict[U] = annihilator(homology(K, k)[1])
  end
  return IdealSheaf(X, ideal_dict; check=false)
end

function common_refinement(list::Vector{<:Covering}, def_cov::Covering)
  isempty(list) && error("list of coverings must not be empty")

  if length(list) == 1
    result = first(list)
    return result, [identity_map(result)]
  end
  patch_list = AbsAffineScheme[]
  anc_list = AbsAffineScheme[]
  to_U_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()
  to_V_dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}()

  if length(list) == 2
    for U in patches(list[1])
      match_found = false
      for V in patches(list[2])
        success, W = _have_common_ancestor(U, V)
        !success && continue
        match_found = true
        push!(anc_list, W)
        #inc_U = _flatten_open_subscheme(U, W)
        #inc_V = _flatten_open_subscheme(V, W)
        inc_U, h_U = _find_chart(U, W)
        inc_U = PrincipalOpenEmbedding(inc_U, h_U; check=false)
        inc_V, h_V = _find_chart(V, W)
        inc_V = PrincipalOpenEmbedding(inc_V, h_V; check=false)

        UV, to_U, to_V = fiber_product(inc_U, inc_V) 
        push!(patch_list, UV)
        to_U_dict[UV] = to_U
        to_V_dict[UV] = to_V
      end
      !match_found && error("no common ancestor found for $U and $V")
    end
    #anc_cov = Covering(anc_list)
    #inherit_glueings!(anc_cov, def_cov)
    result = Covering(patch_list)
    inherit_glueings!(result, def_cov)

    tot_inc1 = CoveringMorphism(result, list[1], to_U_dict; check=false)
    tot_inc2 = CoveringMorphism(result, list[2], to_V_dict; check=false)
    return result, [tot_inc1, tot_inc2]
  end

  # More than two entries
  n = length(list)
  k = div(n, 2)
  res1, inc1 = common_refinement(list[1:k], def_cov)
  res2, inc2 = common_refinement(list[k+1:end], def_cov) 

  result, inc_tot = common_refinement([res1, res2], def_cov)
  return result, vcat([compose(inc_tot[1], inc1[k]) for k in 1:length(inc1)], 
                      [compose(inc_tot[2], inc2[k]) for k in 1:length(inc2)]
                     )
end

function strict_transform(bl::AbsSimpleBlowdownMorphism, inc::ClosedEmbedding)
  B = codomain(bl)
  @assert length(affine_charts(B)) == 1 && first(affine_charts(B)) === codomain(inc)
  inc_cov = CoveredClosedEmbedding(inc, codomain=B)
  return strict_transform(bl, inc_cov)
end

is_graded(R::Ring) = false

########################################################################
# Implement the interface specified in 
# experimental/Schemes/BlowupMorphism.jl
########################################################################

function exceptional_divisor(phi::BlowUpSequence)
  #TODO: Check whether the full resolution has already been computed?
  if !isdefined(phi, :exceptional_divisor)
    phi.exceptional_divisor = sum(phi.ex_div; init=CartierDivisor(domain(phi), ZZ))
  end
  return phi.exceptional_divisor
end

function exceptional_locus(phi::BlowUpSequence)
  if !isdefined(phi, :exceptional_locus)
    phi.exceptional_locus = weil_divisor(exceptional_divisor(phi))
  end
  return phi.exceptional_locus
end

# The following two methods will be ambiguous in general
# so we need to repeat the procedure for the specific types 
# of the second argument.
function strict_transform(phi::BlowUpSequence, a::Any)
  for psi in morhpisms(phi)
    a = strict_transform(psi, a)
  end
  return a
end

function total_transform(phi::BlowUpSequence, a::Any)
  for psi in morphisms(phi)
    a = total_transform(psi, a)
  end
  return a
end

function strict_transform(phi::BlowUpSequence, a::AbsIdealSheaf)
  for psi in morphisms(phi)
    a = strict_transform(psi, a)
  end
  return a
end

function total_transform(phi::BlowUpSequence, a::AbsIdealSheaf)
  for psi in morphisms(phi)
    a = total_transform(psi, a)
  end
  return a
end

function center(phi::BlowUpSequence)
  return center(last_map(phi))
end

