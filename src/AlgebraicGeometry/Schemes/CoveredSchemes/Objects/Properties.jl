
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

function is_smooth_parallel2(X::CoveredScheme)
  channels = Oscar.params_channels(Ring)
  n = length(channels)
  wids = workers()
  for (i, U) in enumerate(affine_charts(X))
    _put_params_rec(channels, OO(U))
  end
  cov = default_covering(X)
  data = [(U, ideal(OO(U), decomposition_info(cov)[U])) for U in affine_charts(X)]
  function compute(dat)
    return _is_smooth(dat[1]; focus=dat[2])
  end
  result = pmap(compute, data)
  return all(isone(x) for x in result)
end

function _is_smooth(U::AbsAffineScheme{<:Field, RT};
    focus::Ideal
  ) where {RT <: MPolyRing}
  return true
end

function _is_smooth(U::AbsAffineScheme{<:Field, RT};
    focus::Ideal
  ) where {RT <: MPolyLocRing}
  return true
end

function _is_smooth(U::AbsAffineScheme{<:Field, RT};
    focus::Ideal = ideal(OO(U), elem_type(OO(U))[])
  ) where {RT <: Union{MPolyQuoRing, MPolyQuoLocRing}}
  if is_known_to_be_irreducible(U) || true
    g = gens(modulus(OO(U)))
    Q, pr = quo(OO(U), focus)
    A = map_entries(Q, jacobian_matrix(g))
    return has_constant_rank(A)
  end
  return error("not implemented")
end

function has_constant_rank(
    A::MatrixElem;
    focus::Union{Ideal, Nothing} = nothing
  )
  @show size(A)
  if focus !== nothing
    Q, pr = quo(base_ring(A), focus)
    AA = map_entries(pr, A)
    return has_constant_rank(AA)
  end
  R = base_ring(A)
  bound = 5
  if nrows(A) <= bound || ncols(A) <= bound
    for r in 1:3
      @show r
      I = ideal(R, minors(A, r))
      is_one(I) && continue
      is_zero(I) && return true
      return false
    end
  end
      
  is_zero(A) && return true
  m = nrows(A)
  n = ncols(A)
  entry_list = elem_type(R)[A[i, j] for i in 1:m for j in 1:n]
  I = ideal(R, entry_list)
  is_one(I) || return false
  c = coordinates(one(R), I)
  ind = _non_zero_indices(c)
  @show length(ind)
  ind_pairs = [(div(k-1, n)+1, mod(k-1, n)+1) for k in ind]
  foc_eqns = elem_type(R)[]
  for (i, j) in ind_pairs
    @show "$(length(foc_eqns)+1)-th entry at ($i, $j) out of $(length(ind_pairs))"
    U = powers_of_element(lifted_numerator(A[i, j]))
    focus = ideal(R, foc_eqns)
    Q, pr = quo(R, focus)
    L, loc_map = localization(Q, U)
    tmp = map_entries(pr, A)
    LA = map_entries(loc_map, tmp)
    multiply_row!(LA, inv(LA[i, j]), i)
    for k in 1:m
      k == i && continue
      add_row!(LA, -LA[k, j], i, k)
    end
    sub = hcat(vcat(LA[1:i-1, 1:j-1], LA[i+1:m, 1:j-1]),
               vcat(LA[1:i-1, j+1:n], LA[i+1:m, j+1:n]))
    @time res = has_constant_rank(sub) 
    res || return false
    push!(foc_eqns, A[i, j])
  end
  return true
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
