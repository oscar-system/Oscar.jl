# An abstract type from which all concrete tasks should be derived. 
abstract type ParallelTask end 

# A generic implementation to extract all parent-like objects appearing 
# in the fields of a concrete instance of a `ParallelTask`. Knowing these 
# beforehand is necessary to be able to create all required parents on 
# the nodes up front. 
function type_params(pt::T) where T <: ParallelTask
  return TypeParams(
    T,
    map(n -> Symbol(n) => type_params(getfield(pt, n)), fieldnames(T))...)
end

function save_type_params(s::SerializerState,
                          tp::TypeParams{T, <: Tuple{Vararg{Pair}}}) where T <: ParallelTask
  save_data_dict(s) do
    save_object(s, encode_type(type(tp)), :name)
    save_data_dict(s, :params) do
      for (k, param_tp) in params(tp)
        save_type_params(s, param_tp, k)
      end
    end
  end
end

# A generic method to create all parents on the node which are required for 
# sending the data in the task's fields. 
function load_type_params(s::DeserializerState, T::Type{<: ParallelTask})
  # If there are no `:params`, quit.
  !haskey(s, :params) && return T, nothing
  params_dict = Dict{Symbol, Any}()
  fields = Type[]
  # Go through the fields of the tasks and load their type parameters.
  load_node(s, :params) do _
    # Note that for this loop we actually need the *concrete* types of the 
    # fields of `T`. That's why they need to be clear from the beginning. 
    for (n,t) in zip(fieldnames(T), fieldtypes(T))
      n == :__attrs && continue # Don't bother with attributes.
      load_node(s, Symbol(n)) do _
        U = decode_type(s)
        field_type, params = load_type_params(s, U)
        push!(fields, field_type)
        params_dict[Symbol(n)] = params
      end
    end
  end

  return T{fields...}, params_dict
end


### Generic methods to send parent-like objects up front. 

# Initial method for sending the parent-like objects in a task. 
function put_type_params(channel::RemoteChannel, obj::T) where T <: ParallelTask
  # `type_params(obj)` returns a pair `(T, params::Dict)` where `T` is the 
  # type itself and `params` is a dictionary with the output of `type_params` 
  # for the fields of `T`.
  for param_tp in params(type_params(obj))
    # For every field go into recursion.
    put_type_params(channel, param_tp.second)
  end
end

# Method for end of recursion.
# TODO: there is probably a way to stop the code from getting here that will speed up the code
function put_type_params(channel::RemoteChannel, ::TypeParams{T, Nothing}) where T
  return
end

function put_type_params(channel::RemoteChannel,
                         ::TypeParams{T, Tuple{Vararg{Pair{Symbol, Nothing}}}}) where T
  return
end

# Recursive call. Send all subsequent parents on which this object is 
# based and finally the object itself, if applicable. 
function put_type_params(channel::RemoteChannel, tp::TypeParams)
  # only  types that use ids need to be sent to the other processes
  put_type_params(channel, params(tp))
end

function put_type_params(channel::RemoteChannel, tps::Tuple{Vararg{Pair}})
  for tp in tps
    put_type_params(channel, tp.second)
  end
end

function put_type_params(channel::RemoteChannel, obj::T) where T
  # only  types that use ids need to be sent to the other processes
  if serialize_with_id(T)
    put_type_params(channel, type_params(obj))
    put!(channel, obj)
  else
    put_type_params(channel, type_params(obj))
  end
end

#function put_type_params(channel::RemoteChannel, v::Vector)
#  for tp in v
#    put_type_params(channel, params(tp))
#  end
#end

### Generic methods to (de-)serialize a concrete task. 
function save_object(s::SerializerState, obj::T) where T <: ParallelTask
  save_data_dict(s) do
    for n in fieldnames(T)
      n == :__attrs && continue
      save_object(s, getfield(obj, n), Symbol(n))
    end
  end
end

function load_object(s::DeserializerState, ::Type{T}, params::Dict) where T <: ParallelTask
  fields = []
  for n in fieldnames(T)
    n == :__attrs && continue
    push!(fields, load_object(s, fieldtype(T, n), params[n], Symbol(n)))
  end
  return T(fields...)
end

# The method of `_compute` for the concrete task specifies what to do 
# on the respective worker. The data will be extracted from the task's 
# fields. The return value must be of the form `(success::Bool, result)` 
# where `success` is to indicate whether we have an affirmative result 
# in any reasonable sense. It is used to decide whether `all` or `any` 
# of a given list of tasks to be done in parallel is achieved. 
#
# The second value `result` can be any reasonable result of the computation 
# which is to be returned to the main node. Note that you must not create 
# new parents on the worker which are required for the contents of `return`, 
# i.e. they need to use the parent-like objects sent from the main node. 
function _compute(::T) where T <: ParallelTask
  error("please implement the function _compute for the type $T")
end

########################################################################
# Generic implementations for deployment of tasks
########################################################################
@doc raw"""
    parallel_all(
        task_list::Vector{TaskType};
        workers::Vector{Int}=Oscar.workers(),
      ) where {TaskType <: ParallelTask}

Given a list `tasklist` of `ParallelTask`s and a pool of workers, deploy them and wait for 
their results. In case all computations were successful, return `(true, res_list)` where 
`res_list` is the `Vector` of results of the respective tasks in the same order. 
If any of the computations was not successful, return `(false, res_list)`.

The user can specify a list of worker ids to be used for deployment via the kwarg `workers`.
"""
function parallel_all(
    task_list::Vector{TaskType};
    workers::Vector{Int}=Oscar.workers(), # Specify which workers to use
    channel_size::Int=32
  ) where {TaskType <: ParallelTask} # TaskType is the type of the task to be deployed.
  n = length(task_list)
  w = length(workers)
  is_zero(w) && !isempty(task_list) && error("zero workers available for non-trivial task; aborting")
  fut_vec = _deploy_work(task_list, workers; channel_size)
  @sync fut_vec
  results = [fetch(fut) for (fut, _) in fut_vec]
  return all(success for (success, _) in results), [result for (_, result) in results]
end

@doc raw"""
    parallel_any(
        task_list::Vector{T};
        workers::Vector{Int}=Oscar.workers(),
        kill_workers::Bool=false
      ) where {T <: ParallelTask}

Given a list `tasklist` of `ParallelTask`s and a pool of workers, deploy them and wait for 
the first affirmative result to come back. In that case return `(true, k, result)` where 
`k` is the number of the successful task `task_list` and `result` its result. 
If all tasks return `(false, _)`, this function returns `(false, 0, nothing)`.

The user can specify a list of worker ids to be used for deployment via the kwarg `workers`.
When `kill_workers` is set to `true`, the workers are killed after call to this function.
"""
function parallel_any(
    task_list::Vector{T};
    workers::Vector{Int}=Oscar.workers(), # Specify which workers to use
    kill_workers::Bool=false,
    channel_size::Int=32
  ) where {T <: ParallelTask} # T is the type of the task to be deployed.
  n = length(task_list)
  w = length(workers)
  is_zero(w) && !isempty(task_list) && error("zero workers available for non-trivial task; aborting")
  fut_vec = _deploy_work(task_list, workers; channel_size)
  while true
    all_failed = true
    for (k, (fut, wid)) in enumerate(fut_vec)
      if !isready(fut)
        all_failed = false # We don't know yet
        continue
      end
      success, result = fetch(fut)
      if success
        # kill the workers if asked for
        kill_workers && map(rmprocs, workers)
        return success, k, result
      end
    end
    all_failed && return false, 0, nothing
  end
end

# Internal method to send tasks for computation to a pool of workers. 
# Returns a `Vector` of `Tuple`s `(fut, wid)` of the `Future`s and the 
# id of the worker where the task has been sent to.
function _deploy_work(
    task_list::Vector{TaskType},
    workers::Vector{Int};
    channel_size::Int=32
  ) where {TaskType}
  w = length(workers)
  individual_channels = Dict{Int, RemoteChannel}(i => RemoteChannel(()->Channel{Any}(channel_size), i) for i in workers)
  assigned_workers = IdDict{TaskType, Int}()
  fut_vec = Tuple{Future, Int}[]
  for (i, task) in enumerate(task_list)
    wid = workers[mod(i, w) + 1]
    channel = individual_channels[wid]
    put_type_params(channel, task)
    assigned_workers[task] = wid
    fut = remotecall(_compute, wid, task)
    push!(fut_vec, (fut, wid))
  end
  return fut_vec
end

