using Distributed: RemoteChannel, Future, remotecall, @everywhere, WorkerPool, AbstractWorkerPool, addprocs, rmprocs, remotecall_eval, nworkers
import Distributed: remotecall, workers, remotecall_fetch

import .Serialization: put_type_params

@doc raw"""
     OscarWorkerPool

An `OscarWorkerPool` is used to handle a pool of separate worker processes
running Oscar. The messaging between processes in an `OscarWorkerPool`
uses the Oscar serialization which is needed to implement parallel methods.
The julia [Distributed](https://docs.julialang.org/en/v1/stdlib/Distributed/)
functions [`remotecall`](https://docs.julialang.org/en/v1/stdlib/Distributed/#Distributed.remotecall-Tuple{Any,%20AbstractWorkerPool,%20Vararg{Any}}),
[remotecall_fetch](https://docs.julialang.org/en/v1/stdlib/Distributed/#Distributed.remotecall_fetch-Tuple{Any,%20AbstractWorkerPool,%20Vararg{Any}})
and [`pmap`](https://docs.julialang.org/en/v1/stdlib/Distributed/#Distributed.pmap) will work with an `OscarWorkerPool` which is a subtype of
[`AbstractWorkerPool`](https://docs.julialang.org/en/v1/stdlib/Distributed/#Distributed.AbstractWorkerPool).
"""
mutable struct OscarWorkerPool <: AbstractWorkerPool
  wp::WorkerPool # the plain worker pool
  channel::Channel{Int} # required for the `AbstractWorkerPool` interface
  workers::Set{Int}     # same
  wids::Vector{Int}     # a list of the ids of all workers which are associated to this pool
  oscar_channels::Dict{Int, <:RemoteChannel} # channels for sending `type_params` to the workers

  function OscarWorkerPool(n::Int)
    wids = addprocs(n)
    wp = WorkerPool(wids)
    # @everywhere can only be used on top-level, so have to do `remotecall_eval` here.
    remotecall_eval(Main, wids, :(using Oscar))

    return new(wp, wp.channel, wp.workers, wids, Dict{Int, RemoteChannel}())
  end
end

@doc raw"""
     oscar_worker_pool(n::Int)
     oscar_worker_pool(f::Function, n::Int)

Create an `OscarWorkerPool` with `n` separate processes running Oscar.
There is also the option to use an `OscarWorkerPool` within a context,
such that closing down the processes happens automatically.

# Example
The following code will start up 3 processes with Oscar,
run a parallel computation over each element in an array and
then shutdown the processes.
```
results = oscar_worker_pool(3) do wp
  Qxy, (x, y) = QQ[:x, :y]
  pmap(z -> z^2, wp, [x^2, x*y, x*y^2])
end
```
"""
oscar_worker_pool(n::Int) = OscarWorkerPool(n)

function oscar_worker_pool(f::Function, n::Int)
  wp = OscarWorkerPool(n)
  local results
  try
    results = f(wp)
  catch e
    rethrow(e)
  finally
    close!(wp)    
  end
  return results
end

### The following implements the `AbstractWorkerPool` interface

# Add a new worker to the pool
function push!(wp::OscarWorkerPool, id::Int)
  # Make sure the node is running Oscar
  remotecall_eval(Main, id, :(using Oscar))
  push!(wp.wids, id) # update the list of associated workers
  return push!(wp.wp, a)
end

# Take a worker from the pool; this marks it as being busy.
function Base.take!(wp::OscarWorkerPool)
  return take!(wp.wp)
end

# Put a formerly busy worker back into the pool.
function Base.put!(wp::OscarWorkerPool, a::Int)
  return put!(wp.wp, a)
end

# Get the number of all workers associated to the pool.
function length(wp::OscarWorkerPool)
  return length(wp.wp)
end

# Return whether any worker is available for work.
function Base.isready(wp::OscarWorkerPool)
  return isready(wp.wp)
end

### end of implementation interface `AbstractWorkerPool`

### extra functionality

workers(wp::OscarWorkerPool) = wp.wids

function get_channel(wp::OscarWorkerPool, id::Int; channel_size::Int=1024)
  chnl = get!(wp.oscar_channels, id) do
    RemoteChannel(()->Channel{Any}(channel_size), id)
  end
end

close!(wp::OscarWorkerPool) = map(rmprocs, workers(wp))

# extend functionality so that `pmap` works with Oscar stuff

function put_type_params(wp::OscarWorkerPool, a::Any)
  for id in wp.wids
    put_type_params(get_channel(wp, id), a)
  end
end

### distributed computation with centralized data management

@doc raw"""
    compute_distributed!(ctx, wp::OscarWorkerPool; wait_period=0.1)

Given an `OscarWorkerPool` `wp` and a context object `ctx`, distribute tasks
coming out of `ctx` to be carried out on the workers (based on availability) and
process the result.

For this to work, the following methods must be implemented for `ctx`:
  - `pop_task!(ctx)` to either return a `Tuple` `(task_id::Int, func::Any, args)`, in which case `remotecall(func, wp, args)` will be called to deploy the task to the workers, or return `nothing` to indicate that all tasks in `ctx` have been exhausted, or return an instance of `WaitForResults` to indicate that some information on results of tasks which have already been given out is needed to proceed with giving out new tasks.
  - `process_result!(ctx, task_id::Int, res)` to process the result `res` of the computation of the task with id `task_id`.

Note: The programmers themselves are responsible to deliver `task_id`s which are unique for the respective tasks!

Deploying tasks to the workers continues until `pop_task!(ctx)` returns `nothing`.
Computation continues until all deployed tasks have been treated with `process_result!`.
The return value is the `ctx` in its current state after computation.
"""
function compute_distributed!(ctx, wp::OscarWorkerPool; wait_period=0.1)
  !isready(wp) && error("pool is not ready")
  is_zero(nworkers(wp)) && error("can not do computations on empty pool")
  n = nworkers(wp)
  futures = Dict{Int, Future}() # task_id => future
  tasks_exhausted = false
  while true
    while isready(wp)
      next_task = pop_task!(ctx)
      if isnothing(next_task)
        tasks_exhausted = true
        break
      elseif next_task isa WaitForResults
        isempty(futures) && error("context object indicates that it is waiting for results, but no `Future`s are present for retrieval")
        break
      end
      task_id, func, args = next_task
      futures[task_id] = remotecall(func, wp, args)
    end

    # gather results
    for (task_id, fut) in futures
      isready(fut) || continue
      process_result!(ctx, task_id, fetch(fut))
      delete!(futures, task_id)
    end

    # if we're done, we're done
    tasks_exhausted && isempty(futures) && break

    # prevent spin-locking
    sleep(wait_period)
  end
  return ctx
end

@doc raw"""
    WaitForResults

A dummy type so that `pop_task!` can return `WaitForResults()` to indicate
that no new task can be deployed before retrieving results of those already
given out.
"""
struct WaitForResults end

@doc raw"""
    pop_task!(ctx::CtxType) where {CtxType}

Internal function for `compute_distributed!`; to be overwritten with methods for specific `CtxType`s.

Whatever type of object the user would like to use as context object `ctx` for `compute_distributed`
should implement a method of this function to do the following.
  - If there is another task to be deployed to a worker, return a triple `(task_id, func, args)` consisting of a function `func`, a unique `task_id::Int`, and arguments `arg` so that `func(arg)` is called on some worker.
  - If no new task can be given out with the current state of the context object `ctx` and we first need to wait for processing of some results of tasks given out already, return `WaitForResults()`.
  - If all tasks in `ctx` have been exhausted, return `nothing`.
"""
pop_task!(ctx) = nothing

@doc raw"""
    process_result!(ctx::CtxType, id::Int, res) where {CtxType}

Internal function for `compute_distributed!`; to be overwritten with methods for specific `CtxType`s.

For a task `(task_id, func, args)` returned by a call to `pop_task!(ctx)`, the result
of the call `func(args)` is delivered as `res` for processing on the main process.
"""
process_result!(ctx, id::Int, res) = ctx

#=
# Custom versions of `remotecall` and `remotecall_fetch` for `OscarWorkerPool`s.
# The idea is to use this special type of worker pool to send the `type_params`
# of the arguments to the workers up front. Then hopefully, whenever an
# `OscarWorkerPool` is being used, loads of other functionality like `pmap`
# will run out of the box.
=#
function remotecall(f::Any, wp::OscarWorkerPool, args...; kwargs...)
  wid = take!(wp)
  for a in args
    put_type_params(get_channel(wp, wid), a)
  end
  for a in kwargs
    put_type_params(get_channel(wp, wid), a)
  end

  # Copied from Distributed.jl/src/workerpool.jl.
  # This puts the worker back to the pool once the future is ready
  # and we do not have to worry about this ourselves.
  local fut
  try
    fut = remotecall(f, wid, args...; kwargs...)
  catch
    put!(wp, wid)
    rethrow()
  end
  t = Threads.@spawn Threads.threadpool() try
    wait(fut)
  catch
  finally
    put!(wp, wid)
  end
  errormonitor(t)

  return fut
end

function remotecall_fetch(f::Any, wp::OscarWorkerPool, args...; kwargs...)
  wid = take!(wp)
  for a in args
    put_type_params(get_channel(wp, wid), a)
  end
  for a in kwargs
    put_type_params(get_channel(wp, wid), a)
  end
  result = remotecall_fetch(f, wid, args...; kwargs...)
  put!(wp, wid)
  return result
end

export oscar_worker_pool
