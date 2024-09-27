
########################################################################
# Properties of AbsCoveredSchemes                                      #
#    by referral to the patches                                        #
########################################################################

########################################################################
# (1) Check and store, whether a covered scheme is empty               #
########################################################################
@doc raw"""
    is_empty(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is empty.

"""
is_empty(X::AbsCoveredScheme) = is_empty(underlying_scheme(X))
@attr Bool function is_empty(X::CoveredScheme)
  if !isdefined(X, :coverings)
    return true
  end
  return all(isempty, all_patches(default_covering(X)))
end

########################################################################
# (2) Check and store, whether a covered scheme is smooth              #
########################################################################
@doc raw"""
    is_smooth(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is smooth.
"""
is_smooth(X::AbsCoveredScheme) = is_smooth(underlying_scheme(X))

@attr Bool function is_smooth(X::CoveredScheme)
  if !isdefined(X, :coverings)
    return true
  end
  # TODO: Can this be optimized? We only need to check the rank of
  # the Jacobian matrices where we haven't checked in another chart before.
  # This is a tradeoff between cost of Jacobian criterion and cost of
  # optimizing the use of the covering
  return all(is_smooth, affine_charts(X))
end

### Recursive function to broadcast parent-like objects to all other workers
function _put_params_rec(channels::Vector{RemoteChannel{Channel{T}}}, a::Any) where {T}
  !serialize_with_id(a) && return # check whether we actually need to send something; if not, quit.
  error("recursive broadcasting is not implemented for objects of type $(typeof(a))")
end

function _put_params_rec(channels::Vector{RemoteChannel{Channel{T}}}, a::MPolyRing) where {T}
  _put_params_rec(channels, coefficient_ring(a))
  put_params(channels, a)
  @hassert :Parallelization 2 _check_other_sides(a)
end

function _put_params_rec(channels::Vector{RemoteChannel{Channel{T}}}, a::MPolyQuoRing) where {T}
  _put_params_rec(channels, base_ring(a))
  put_params(channels, a)
  @hassert :Parallelization 2 _check_other_sides(a)
end

function _put_params_rec(channels::Vector{RemoteChannel{Channel{T}}}, a::MPolyLocRing) where {T}
  _put_params_rec(channels, base_ring(a))
  put_params(channels, inverted_set(a))
  put_params(channels, a)
  @hassert :Parallelization 2 _check_other_sides(a)
end

function _put_params_rec(channels::Vector{RemoteChannel{Channel{T}}}, a::MPolyQuoLocRing) where {T}
  _put_params_rec(channels, base_ring(a))
  _put_params_rec(channels, localized_ring(a))
  _put_params_rec(channels, underlying_quotient(a))
  put_params(channels, a)
  @hassert :Parallelization 2 _check_other_sides(a)
end

function _check_other_sides(R::Any)
  ref = get(global_serializer_state.obj_to_id, R, nothing)
  ref === nothing && return false
  for wid in workers()
    remotecall_fetch((a)->haskey(Oscar.global_serializer_state.id_to_obj, Oscar.UUID(a)), wid, string(ref)) || return false
  end
  return true
end

### Functions to make parent-like objects known to particular 
# other workers through particular channels
function _put_and_take(channel::RemoteChannel, wid::Int, a::Any)
  @show a
  @show serialize_with_id(a)
  !serialize_with_id(a) && return # check whether we actually need to send something.
  @show length(keys(global_serializer_state.obj_to_id))
  put!(channel, a)
  @show length(keys(global_serializer_state.obj_to_id))
  function take_param(ch::RemoteChannel)
    params = take!(ch)
  end
  remotecall_wait(take_param, wid, channel)
  @show length(keys(global_serializer_state.obj_to_id))
  _check_other_side(a, wid)
end

function _put_ring(channel::RemoteChannel{T}, wid::Int, R::MPolyRing) where {T<:Channel{<:Ring}}
  @show "sending polynomial ring"
  _put_and_take(channel, wid, coefficient_ring(R))
  _put_and_take(channel, wid, R)
  @show "done sending polynomial ring"
end

function _put_ring(channel::RemoteChannel{T}, wid::Int, R::MPolyQuoRing) where {T<:Channel{<:Ring}}
  _put_ring(channel, wid, base_ring(R))
  _put_and_take(channel, wid, R)
end

function _put_ring(channel::RemoteChannel{T}, wid::Int, R::MPolyLocRing) where {T<:Channel{<:Ring}}
  _put_ring(channel, wid, base_ring(R))
  _put_and_take(channel, wid, R)
end

function _put_ring(channel::RemoteChannel{T}, wid::Int, R::MPolyQuoLocRing) where {T<:Channel{<:Ring}}
  _put_ring(channel, wid, base_ring(R))
  _put_ring(channel, wid, underlying_quotient(R))
  _put_ring(channel, wid, localized_ring(R))
  _put_and_take(channel, wid, R)
end


function _check_other_side(R::Any, wid::Int)
  ref = get(global_serializer_state.obj_to_id, R, nothing)
  @show typeof(R)
  @assert ref !== nothing
  @assert remotecall_fetch((a)->haskey(Oscar.global_serializer_state.id_to_obj, Oscar.UUID(a)), wid, string(ref))
end

function _compute_(dat)
  return _is_smooth(dat[1]; focus=dat[2])
end

function is_smooth_parallel2(X::CoveredScheme)
  channels = Oscar.params_channels(Ring)
  n = length(channels)
  wids = workers()
  for (i, U) in enumerate(affine_charts(X))
    _put_params_rec(channels, OO(U))
  end
  cov = default_covering(X)
  data = [(U, ideal(OO(U), decomposition_info(cov)[U])) for U in affine_charts(X)]
  result = pmap(_compute_, data)
  return all(isone(x) for x in result)
end

function _is_smooth_parallel(X::AbsAffineScheme)
  return _is_smooth_parallel(underlying_scheme(X))
end

function _is_smooth(U::AbsAffineScheme{<:Field, RT};
    focus::Ideal = ideal(OO(U), elem_type(OO(U))[]),
    jacobian_cut_off::Int=5
  ) where {RT <: MPolyRing}
  return true
end

function _is_smooth(U::AbsAffineScheme{<:Field, RT};
    focus::Ideal = ideal(OO(U), elem_type(OO(U))[]),
    jacobian_cut_off::Int=5
  ) where {RT <: MPolyLocRing}
  return true
end

function _is_smooth(U::AbsAffineScheme{<:Field, RT};
    focus::Ideal = ideal(OO(U), elem_type(OO(U))[]),
    jacobian_cut_off::Int=5
  ) where {RT <: Union{MPolyQuoRing, MPolyQuoLocRing}}
  g = gens(modulus(OO(U)))
  Q, pr = quo(OO(U), focus)
  A = map_entries(Q, transpose(jacobian_matrix(g)))
  return has_locally_constant_corank(A; jacobian_cut_off, upper_bound=dim(U))[1]
end

function _is_smooth_parallel(U::AffineScheme{<:Field, RT};
    focus::Ideal = ideal(OO(U), elem_type(OO(U))[]),
    jacobian_cut_off::Int=5
  ) where {RT <: Union{MPolyQuoRing, MPolyQuoLocRing}}
  g = gens(modulus(OO(U)))
  Q, pr = quo(OO(U), focus)
  A = map_entries(Q, transpose(jacobian_matrix(g)))
  return has_locally_constant_corank_parallel(A; jacobian_cut_off, upper_bound=dim(U))[1]
end

_complexity(a::RingElem) = 0
_complexity(a::MPolyRingElem) = length(a)
_complexity(a::MPolyQuoRingElem) = length(lift(a))
_complexity(a::MPolyLocRingElem) = length(numerator(a)) + length(denominator(a))
_complexity(a::MPolyQuoLocRingElem) = length(lifted_numerator(a)) + length(lifted_denominator(a))


function has_locally_constant_corank_parallel(
    A::MatrixElem;
    upper_bound::Union{Int, Nothing}=nothing,
    jacobian_cut_off::Int=5
  )
  is_zero(A) && return true, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]

  R = base_ring(A)
  if is_zero(ngens(R))
    AA = map_entries(x->constant_coefficient(lifted_numerator(x))*inv(constant_coefficient(lifted_denominator(x))), A)
    return true, [(ideal(R, elem_type(R)[]), ncols(A) - rank(AA))]
  end

  m = nrows(A)
  n = ncols(A)
  entry_list = elem_type(R)[A[i, j] for i in 1:m for j in 1:n]
  I = ideal(R, entry_list)
  if !is_one(I) 
    # It might still be the case that we have different ranks on different components; see above
    J = quotient(ideal(R, elem_type(R)[]), I)
    is_one(I + J) || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
    Q1, pr1 = quo(R, I)
    res1, list1 = has_locally_constant_corank_parallel(map_entries(pr1, A); upper_bound, jacobian_cut_off)
    res1 || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]

    Q2, pr2 = quo(R, J)
    res2, list2 = has_locally_constant_corank_parallel(map_entries(pr2, A); upper_bound, jacobian_cut_off)
    return res2, vcat([(ideal(R, [R(g) for g in gens(saturated_ideal(I1))]) + I, k) for (I1, k) in list1],
                      [(ideal(R, [R(g) for g in gens(saturated_ideal(I2))]) + J, k) for (I2, k) in list2])
  end
  c = coordinates(one(R), I)
  ind = _non_zero_indices(c)
  ind_with_comp = [(k, _complexity(c[ind[k]])) for k in 1:length(ind)]
  sort!(ind_with_comp; by=p->p[2])
  ind = [ind[k] for k in first.(ind_with_comp)]
  #@show length(ind)
  ind_pairs = [(div(k-1, n)+1, mod(k-1, n)+1) for k in ind]
  foc_eqns = elem_type(R)[]
  full_list = Vector{Tuple{Ideal, Int}}()
  data = Tuple[]
  maps = Tuple[]
  simplification = false
  for (i, j) in ind_pairs
    simplification = false
    #@show "$(length(foc_eqns)+1)-th entry at ($i, $j) out of $(length(ind_pairs))"
    U = powers_of_element(lifted_numerator(A[i, j]))
    focus = ideal(R, foc_eqns)
    focus = ideal(R, small_generating_set(focus)) # keep the data small
    #@show is_one(focus)
    is_one(focus) && break
    Q, pr = quo(R, focus)
    L, loc_map = localization(Q, U)
    L_simp, simp_map, simp_map_inv = simplify(L)
    if ngens(L_simp) < ngens(L) 
      simplification = true
      push!(maps, (simp_map, simp_map_inv))
    else
      push!(maps, (nothing, nothing))
    end
    tmp = map_entries(pr, A)
    help_map = x->L(lifted_numerator(x), lifted_denominator(x); check=false)
    LA = map_entries(help_map, tmp)
    simplification && (LA = map_entries(simp_map, tmp2))
    LA = map_entries(_simplify, LA) # keep the data to be sent small
    u = LA[i, j]
    #multiply_row!(LA, inv(LA[i, j]), i)
    for k in 1:m
      k == i && continue
      c = -LA[k, j]
      multiply_row!(LA, u, k)
      add_row!(LA, c, i, k)
    end
    sub = hcat(vcat(LA[1:i-1, 1:j-1], LA[i+1:m, 1:j-1]),
               vcat(LA[1:i-1, j+1:n], LA[i+1:m, j+1:n]))
    push!(data, (LA, upper_bound, jacobian_cut_off))
    push!(foc_eqns, A[i, j])
    A[i, j] = zero(R)
  end
  channels = params_channels(Any) # TODO: The union type is necessary here. Why?
  for (LA, _, _) in data
    #@show base_ring(LA)
    _put_params_rec(channels, base_ring(LA))
    put_params(channels, parent(LA))
  end
  results = pmap(_helper_func, data)
  for (i, (success, list)) in enumerate(results)
    !success && return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
    for (J, k) in list
      simp_map_inv = maps[i][2]
      JJ = simp_map_inv !== nothing ? simp_map_inv(J) : J
      J_sat = saturated_ideal(JJ)
      caught = false
      for (K, r) in full_list
        r != k && !is_one(J_sat + K) && return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
        if all(radical_membership(g, K) for g in gens(J_sat))
          caught = true
          break
        end
      end
      !caught && push!(full_list, (J_sat, k))
    end
  end
  return true, [(ideal(R, elem_type(R)[R(g) for g in gens(I)]), k) for (I, k) in full_list]
end

_helper_func(dat) = has_locally_constant_corank(dat[1]; upper_bound=dat[2], jacobian_cut_off=dat[3])
_simplify(a::MPolyRingElem) = a
_simplify(a::MPolyLocRingElem) = a
_simplify(a::MPolyQuoRingElem) = simplify(a)
_simplify(a::MPolyQuoLocRingElem) = parent(a)(lift(simplify(numerator(a))), lifted_denominator(a); check=false)


function has_locally_constant_corank(
    A::MatrixElem;
    upper_bound::Union{Int, Nothing}=nothing,
    jacobian_cut_off::Int=5
  )
  #@show size(A)
  @show "checking for zero matrix..."
  is_zero(A) && return true, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
  @show "done checking for zero matrix"

  R = base_ring(A)
  if is_zero(ngens(R))
    AA = map_entries(x->constant_coefficient(lifted_numerator(x))*inv(constant_coefficient(lifted_denominator(x))), A)
    return true, [(ideal(R, elem_type(R)[]), ncols(A) - rank(AA))]
  end

  if nrows(A) <= jacobian_cut_off || ncols(A) <= jacobian_cut_off
    # use the jacobian criterion
    r0 = upper_bound === nothing ? 1 : ncols(A) - upper_bound
    for r in r0:jacobian_cut_off
      #@show r
      #@show upper_bound
      #@show ncols(A)
      @show "computing minors..."
      I = ideal(R, minors(A, r))
      @show "done computing minors"
      if upper_bound !== nothing && ncols(A) - r + 1 > upper_bound 
        is_one(I) || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
        continue
      end
      #@show is_one(I)
      #@show is_zero(I)
      @show "checking triviality..."
      is_one(I) && continue
      @show "done checking triviality"
      @show "reducing generators..."
      is_zero(I) && return true, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A) - r + 1)]
      @show "done reducing generators..."
      # It could still be the case that we have two disjoint components X and Y 
      # on which A has different ranks. In this case, I will vanish on one of 
      # the components, say X, while 0:I vanishes on Y. 
      J = quotient(ideal(R, elem_type(R)[]), I)
      is_one(I + J) || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]# components are not disjoint in this case
      Q1, pr1 = quo(R, I)
      res1, list1 = has_locally_constant_corank(map_entries(pr1, A); upper_bound, jacobian_cut_off)
      res1 || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]

      Q2, pr2 = quo(R, J)
      res2, list2 = has_locally_constant_corank(map_entries(pr2, A); upper_bound, jacobian_cut_off)
      return res2, vcat([(ideal(R, [preimage(pr1, g) for g in gens(I1)]) + I, k) for (I1, k) in list1],
                        [(ideal(R, [preimage(pr2, g) for g in gens(I2)]) + J, k) for (I2, k) in list2])
    end
  end
      
  m = nrows(A)
  n = ncols(A)
  entry_list = elem_type(R)[A[i, j] for i in 1:m for j in 1:n]
  I = ideal(R, entry_list)
  if !is_one(I) 
    # It might still be the case that we have different ranks on different components; see above
    J = quotient(ideal(R, elem_type(R)[]), I)
    is_one(I + J) || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
    Q1, pr1 = quo(R, I)
    res1, list1 = has_locally_constant_corank(map_entries(pr1, A); upper_bound, jacobian_cut_off)
    res1 || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]

    Q2, pr2 = quo(R, J)
    res2, list2 = has_locally_constant_corank(map_entries(pr2, A); upper_bound, jacobian_cut_off)
    return res2, vcat([(ideal(R, [R(g) for g in gens(saturated_ideal(I1))]) + I, k) for (I1, k) in list1],
                      [(ideal(R, [R(g) for g in gens(saturated_ideal(I2))]) + J, k) for (I2, k) in list2])
  end
  @show "computing coordinates..."
  c = coordinates(one(R), I)
  @show "done computing coordinates"
  ind = _non_zero_indices(c)
  ind_with_comp = [(k, _complexity(c[ind[k]])) for k in 1:length(ind)]
  sort!(ind_with_comp; by=p->p[2])
  ind = [ind[k] for k in first.(ind_with_comp)]
  #@show length(ind)
  ind_pairs = [(div(k-1, n)+1, mod(k-1, n)+1) for k in ind]
  foc_eqns = elem_type(R)[]
  full_list = Vector{Tuple{Ideal, Int}}()
  simplification = false
  for (i, j) in ind_pairs
    simplification = false
    @show "$(length(foc_eqns)+1)-th entry at ($i, $j) out of $(length(ind_pairs))"
    U = powers_of_element(lifted_numerator(A[i, j]))
    focus = ideal(R, foc_eqns)
    is_one(focus) && break
    @show "preparing sub-matrix..."
    Q, pr = quo(R, focus)
    L, loc_map = localization(Q, U)
    L_simp, simp_map, simp_map_inv = simplify(L)
    @show ngens(L_simp) < ngens(L)
    ngens(L_simp) < ngens(L) && (simplification = true)
    tmp = map_entries(pr, A)
    help_map = x->L(lifted_numerator(x), lifted_denominator(x); check=false)
    LA = map_entries(help_map, tmp)
    simplification && (LA = map_entries(simp_map, tmp2))
    multiply_row!(LA, inv(LA[i, j]), i)
    for k in 1:m
      k == i && continue
      add_row!(LA, -LA[k, j], i, k)
    end
    sub = hcat(vcat(LA[1:i-1, 1:j-1], LA[i+1:m, 1:j-1]),
               vcat(LA[1:i-1, j+1:n], LA[i+1:m, j+1:n]))
    @show "done preparing sub-matrix"
    res, list = has_locally_constant_corank(sub; upper_bound, jacobian_cut_off) 
    @show "done with recursive call"
    res || return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
    for (J, k) in list
      JJ = simplification ? simp_map_inv(J) : J
      J_sat = saturated_ideal(JJ)
      caught = false
      for (K, r) in full_list
        r != k && !is_one(J_sat + K) && return false, [(ideal(base_ring(A), elem_type(base_ring(A))[]), ncols(A))]
        if all(radical_membership(g, K) for g in gens(J_sat))
          caught = true
          break
        end
      end
      !caught && push!(full_list, (J_sat, k))
    end
    push!(foc_eqns, A[i, j])
    A[i, j] = zero(R)
  end
  return true, [(ideal(R, elem_type(R)[R(g) for g in small_generating_set(I)]), k) for (I, k) in full_list]
end

_non_zero_indices(c::SRow) = c.pos
_non_zero_indices(c::Vector) = [k for k in 1:length(c) if !is_zero(c[k])]
_non_zero_indices(c::MatrixElem) = [k for k in 1:length(c) if !is_zero(c[1, k])]


function is_known_to_be_irreducible(U::AbsAffineScheme{<:Field, RT}) where {RT <: MPolyRing}
  return true
end

function is_known_to_be_irreducible(U::AbsAffineScheme{<:Field, RT}) where {RT <: MPolyLocRing}
  return true
end

function is_known_to_be_irreducible(
    U::AbsAffineScheme{<:Field, RT}
  ) where {RT <: Union{MPolyQuoRing, MPolyQuoLocRing}}
  has_attribute(U, :is_irreducible) && return is_irreducible(U)
  I = modulus(OO(U))
  if has_attribute(I, :is_prime)
    set_attribute!(U, :is_irreducible=>is_prime(I))
    return is_prime(I)
  end
  return false
end

function is_smooth_parallel(X::CoveredScheme)
  if !isdefined(X, :coverings)
    return true
  end
  # Make sure all required parents are known on the other workers.
  # This step will be automated, eventually, but right now we have 
  # to take care of it.
  channels = Oscar.params_channels(Ring)
  n = length(channels)
  wids = workers()
  @show n
  charts = affine_charts(X)
  data_buckets = Dict{Int, Vector{AbsAffineScheme}}()
  for (i, U) in enumerate(charts)
    @show i
    k = mod(i-1, n)+1
    # pmap chooses freely which worker to use for which entry.
    # Hence, we have no insurance about sending specific rings 
    # to specific workers and, for the moment, we have to send 
    # everything everywhere. 
    for (ch, wid) in zip(channels, wids)
      R = OO(U)
      @show R
      @show typeof(R)
      _put_ring(ch, wid, R)
      ref = get(global_serializer_state.obj_to_id, R, nothing)
      @assert ref !== nothing
      @assert remotecall_fetch((a)->haskey(Oscar.global_serializer_state.id_to_obj, Oscar.UUID(a)), wid, string(ref))
    end
    #=
    @show k
    current_channel = channels[k]
    wid = wids[k]
    _put_ring(current_channel, wid, OO(U))

    # make sure the ring really got to the other side
    R = OO(U)
    ref = get(global_serializer_state.obj_to_id, R, nothing)
    @assert ref !== nothing
    @assert remotecall_fetch((a)->haskey(Oscar.global_serializer_state.id_to_obj, Oscar.UUID(a)), wid, string(ref))
    =#

    bucket = get!(data_buckets, k) do
      AbsAffineScheme[]
    end
    push!(bucket, U)
  end
  
  lst = charts[1:n]
  pmap(one, OO.(lst))


  @show "done with distribution"

  round = 1
  @show data_buckets
  while length(data_buckets) == n && !any(isempty(v) for (_, v) in data_buckets)
    @show "computing round $round"
    round += 1
    data = [pop!(data_buckets[i]) for i in 1:n]

    @show "test whether the rings are over there"
    for (j, U) in enumerate(data)
      @show j
      # make sure the ring really got to the other side
      R = OO(U)
      ref = get(global_serializer_state.obj_to_id, R, nothing)
      @assert ref !== nothing
      @assert remotecall_fetch((a)->haskey(Oscar.global_serializer_state.id_to_obj, Oscar.UUID(a)), wids[j], string(ref))
    end

    @show data
    @show "calling pmap"
    result = pmap(is_smooth, data)
    @show "result: $result"
    any(is_zero(x) for x in result) && return false
  end

  # the last few can be done here
  return all(is_smooth(first(l)) for (_, l) in data_buckets)
end

function _jacobian_criterion(X::CoveredScheme{<:Field})
  if !isdefined(X, :coverings)
    return true
  end

  if !has_decomposition_info(default_covering(X))
    throw(NotImplementedError(:_jacobian_criterion, "only implemented when decomposition info is available"))
  end

  dec_info = decomposition_info(default_covering(X))
  for (V, fs) in dec_info
    R_quotient = OO(V)
    R = R_quotient isa MPolyRing ? R_quotient : base_ring(R_quotient)
    I = saturated_ideal(defining_ideal(V))
    mat = jacobian_matrix(R, gens(I))
    sing_locus = ideal(R, fs) + ideal(R, minors(mat, codim(V)))
    sing_subscheme = subscheme(V, sing_locus)
    if !isempty(sing_subscheme)
      return false
    end
  end
  return true
end

########################################################################
# (3) Check and store, whether a covered scheme is integral            #
#     i.e. irreducible and reduced                                     #
########################################################################
@doc raw"""
    is_integral(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is integral.

"""
@attr Bool function is_integral(X::AbsCoveredScheme)
  !is_empty(X) || return false
  return is_reduced(X) && is_irreducible(X)
end

# auxiliary function for connectedness of gluing graph
#      do not confuse with connectedness of the scheme
# Note: This does not work with gluing_graph, because empty patches
#      need to be ignored without changing the covering
@doc raw"""
    is_connected_gluing(X::AbsCoveredScheme)

Return the boolean value whether the gluing graph of the default
covering of the scheme X is connected.

!!! note
    This function is designed to ignore empty patches, which may arise e.g. upon creation of subschemes of covered schemes,

"""
@attr Bool function is_connected_gluing(X::AbsCoveredScheme)
  return is_connected(pruned_gluing_graph(default_covering(X)))
end

########################################################################
# (4) Check and store, whether a covered scheme is connected           #
########################################################################
@doc raw"""
    is_connected(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is connected.

"""
@attr Bool function is_connected(X::AbsCoveredScheme)
  is_connected_gluing(X) || return false
  # note for future implementation: expensive property
  # 1) do primary decomposition
  # 2) check connectedness of lowest two layers of the intersection lattice
  error("not implemented yet")
end

##############################################################################
# (5) Check and store, whether a scheme is reduced
##############################################################################
@doc raw"""
    is_reduced(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is reduced.

"""
@attr Bool function is_reduced(X::AbsCoveredScheme)
  return all(is_reduced, affine_charts(X))
end

##############################################################################
# (6) Check and store, whether a covered scheme is irreducible
##############################################################################
@doc raw"""
    is_irreducible(X::AbsCoveredScheme)

Return the boolean value whether a covered scheme `X` is irreducible.

"""
@attr Bool function is_irreducible(X::AbsCoveredScheme)
  is_connected_gluing(X) || return false
  !is_empty(X) || return false
  v=findall(!is_empty, affine_charts(X))  ## only check non-empty patches
  return all(is_irreducible(affine_charts(X)[v[i]]) for i in 1:length(v) )
end

# Overwrite computations of determinants as to avoid massive is_zero checks
# and slow arithmetic.
function det(A::MatrixElem{<:MPolyQuoRingElem})
  AA = map_entries(lift, A)
  @assert base_ring(AA) isa MPolyRing
  res = det(AA)
  return base_ring(A)(res; check=false)
end

function det(A::MatrixElem{<:MPolyQuoLocRingElem})
  if all(is_one(lifted_denominator(x)) for x in A)
    AAA = map_entries(lifted_numerator, A)
    res = det(AAA)
    return base_ring(A)(res; check=false)
  end

  AA = map_entries(fraction, A)
  res = det(AA)
  return base_ring(A)(res; check=false)
end

function det(A::MatrixElem{<:MPolyLocRingElem})
  if all(is_one(denominator(x)) for x in A)
    AAA = map_entries(numerator, A)
    res = det(AAA)
    return base_ring(A)(res; check=false)
  end

  AA = map_entries(fraction, A)
  res = det(AA)
  return base_ring(A)(res; check=false)
end


