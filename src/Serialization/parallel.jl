########################################################################
# Simple single-machine parallelization
#
# In this file we explain and implement patterns for parallelization 
# on a single machine with multiple cores. 
#
# This can be used to deploy embarassingly parallel tasks on multiple 
# cores and either wait for all of them to finish, or for the first 
# successful computation on one of the cores. 
########################################################################

# In order to use this infrastructure, you first need to set up a simple 
# record of the data needed to accomplish any concrete instance to a 
# task. 
#
# Say you would like to call `gcd(g...)` where `g` is a list of 
# `RingElem`s, but actually for a whole list `l` of different 
# combinations. Then your record would simply wrap up the arguments 
# you pass to the function, in this case a `Vector` of `RingElem`s.
abstract type ParallelTask end 

function type_params(pt::T) where T <: ParallelTask
  type_params_dict = Dict{Symbol, Any}()
  for n in fieldnames(T)
    if n != :__attrs
      type_params_dict[Symbol(n)] = type_params(getfield(pt, n))
    end
  end
  return typeof(pt), type_params_dict
end

function load_type_params(s::DeserializerState, T::Type{<: ParallelTask})
  !haskey(s, :params) && return T, nothing
  params_dict = Dict{Symbol, Any}()
  fields = DataType[]
  load_node(s, :params) do _
    for (n,t) in zip(fieldnames(T), fieldtypes(T))
      if n!= :__attrs
        load_node(s, Symbol(n)) do _
          U = decode_type(s)
          params = load_type_params(s, U)
          push!(fields, params[1])
          params_dict[Symbol(n)] = params
        end
      end
    end
  end

  return T{fields...}, params_dict
end

function _compute(::T) where T <: ParallelTask
  error("please implement the function _compute for the type $T")
end

function save_object(s::SerializerState, obj::T) where T <: ParallelTask
  save_data_dict(s) do
    for n in fieldnames(T)
      if n != :__attrs
        save_object(s, getfield(obj, n), Symbol(n))
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{T}, params::Dict) where T <: ParallelTask
  fields = []

  for n in fieldnames(T)
    if n!= :__attrs
      push!(fields, load_object(s, params[n]..., Symbol(n)))
    end
  end
  return T(fields...)
end


# The data will need to be passed to the different workers. 
# To allow for this, you need to specify how to serialize your record 
# struct. In many cases, this can be done by the generic serialization 
# implementation. BUT: you still need to implement the following 
# function, which communicates which parent-like object appear in your
# record.

# As said above, the generic (de-)serialization should do, but in case 
# you want something specialized, you can overwrite the methods here. 
#function save_object(s::SerializerState, ds::SampleDataStruct)
#  save_data_dict(s) do
#    save_object(s, ds.elems, :elems)
#  end
#end

#function load_object(s::DeserializerState, ::Type{<:SampleDataStruct}, params::Dict)
#  R = params[:parent]
#  return SampleDataStruct(load_object(s, Vector{elem_type(R)}, R, :elems))
#end

# Finally, implement the function which should be called on an instance 
# of your record (on the workers) to carry out the actual task.
#
# Note that if you want to use `wait_first_parallel`, then the return value must 
# be of the form `(success, result)` where `success` is a `Bool` indicating 
# whether the worker has obtained an affirmative result in any reasonable sense. 
# Execution is stopped only if a return pair with `success==true` is found. 

########################################################################
# Generic implementations for deployment of tasks
########################################################################
function wait_all_parallel(
    task_list::Vector{TaskType};
    workers::Vector{Int}=Oscar.workers(), # Specify which workers to use
  ) where {TaskType} # TaskType is the type of the task to be deployed.
  n = length(task_list)
  w = length(workers)
  is_zero(w) && !isempty(task_list) && error("zero workers available for non-trivial task; aborting")
  fut_vec = _collect_futures(_deploy_work(task_list, workers))
  @sync fut_vec
  return fetch.(fut_vec)
end

function wait_first_parallel(
    task_list::Vector{T};
    workers::Vector{Int}=Oscar.workers(), # Specify which workers to use
    kill_workers::Bool=false
  ) where {T} # T is the type of the task to be deployed.
  n = length(task_list)
  w = length(workers)
  is_zero(w) && !isempty(task_list) && error("zero workers available for non-trivial task; aborting")
  futures = _deploy_work(task_list, workers)
  while true
    for (wid, fut_vec) in futures
      if any(isready, fut_vec)
        k = findfirst(isready, fut_vec)
        fut = fut_vec[k]
        success, result = fetch(fut)
        if success
          # kill the workers if asked for
          kill_workers && map(rmprocs, workers)
          return result
        end
      end
    end
  end
end

function _deploy_work(
    task_list::Vector{TaskType},
    workers::Vector{Int}
  ) where {TaskType}
  w = length(workers)
  individual_channels = Dict{Int, RemoteChannel}(i => RemoteChannel(()->Channel{Any}(32), i) for i in workers)
  assigned_workers = IdDict{TaskType, Int}()
  futures = Dict{Int, Vector{Future}}()
  fut_vec = Future[]
  for (i, task) in enumerate(task_list)
    wid = workers[mod(i, w) + 1]
    channel = individual_channels[wid]
    put!(channel, type_params(task)[2])
    #remotecall(take!, wid, channel)
    assigned_workers[task] = wid
    fut = remotecall(_compute, wid, task)
    push!(get!(futures, wid, Future[]), fut)
    push!(fut_vec, fut)
  end
  return futures
end

function _collect_futures(futures::Dict{Int, Vector{Future}})
  return vcat([futs for (_, futs) in futures]...)
end

