using Distributed: RemoteChannel, Future, remotecall, @everywhere, AbstractWorkerPool, addprocs, rmprocs, remotecall_eval, nworkers, procs
import Distributed: remotecall, workers, remotecall_fetch, pmap, ClusterManager

import .Serialization: put_type_and_params

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

A worker that dies during a computation (it crashed, or was killed by
[`remotecall_with_timeout`](@ref)) is replaced by a new one, if the pool
knows how to start one (see [`oscar_worker_pool`](@ref)).
"""
mutable struct OscarWorkerPool <: AbstractWorkerPool
  # `channel` (free workers) and `workers` (all workers) are the fields the default
  # `AbstractWorkerPool` methods of Distributed work with (`take!`, `put!`, `push!`,
  # `length`, `isready`); in particular `take!` drops dead workers lazily.
  channel::Channel{Int}
  workers::Set{Int}
  oscar_channels::Dict{Int, RemoteChannel} # channels for sending `type_and_params` to the workers
  init_expr::Expr # evaluated on every worker after `using Oscar`
  addworker::Union{Function, Nothing} # () -> id of a new worker; `nothing`: dead workers are not replaced

  function OscarWorkerPool(wids::Vector{Int}, init_expr::Expr, addworker::Union{Function, Nothing})
    wp = new(Channel{Int}(typemax(Int)), Set{Int}(), Dict{Int, RemoteChannel}(), init_expr, addworker)
    # @everywhere can only be used on top-level, so `push!` does `remotecall_eval`.
    asyncmap(wid -> push!(wp, wid), wids)
    return wp
  end
end

function OscarWorkerPool(n::Int, init_expr::Expr=:(); kw...)
  return OscarWorkerPool(addprocs(n; kw...), init_expr, () -> only(addprocs(1; kw...)))
end

function OscarWorkerPool(manager::ClusterManager, init_expr::Expr=:();
                         project=Base.active_project(), addworker=nothing, kw...)
  wids = addprocs(manager; kw...)
  if !isnothing(project)
    asyncmap(wids) do wid
      remotecall_eval(Main, wid, :(using Pkg; Pkg.activate($project); Pkg.instantiate()))
    end
  end
  return OscarWorkerPool(wids, init_expr, addworker)
end

@doc raw"""
     oscar_worker_pool(n::Int, init_expr::Expr=:(); kw...)
     oscar_worker_pool(f::Function, n::Int, init_expr::Expr=:(); kw...)
     oscar_worker_pool(manager::ClusterManager, init_expr::Expr=:(); project=Base.active_project(), addworker=nothing, kw...)
     oscar_worker_pool(f::Function, manager::ClusterManager, init_expr::Expr=:(); project=Base.active_project(), addworker=nothing, kw...)

Create an `OscarWorkerPool` with `n` separate processes running Oscar.
There is also the option to use an `OscarWorkerPool` within a context,
such that closing down the processes happens automatically.
The `init_expr` is evaluated on every worker after loading Oscar, including
workers added to the pool later.

Keyword arguments will get passed to `addprocs` when initializing the workers,
for example use `exeflags` to specify settings for the worker processes.
When passing a `manager`, the `project` argument will determine which project to activate on the machines in the cluster.
The default is to use `Base.active_project()` mimicking the local setting, passing nothing won't activate any project.

A worker that dies is replaced by a new one. With `n` workers the replacement
is started like the others, with a `manager` only if `addworker` is given: a
function without arguments that starts one worker and returns its id, for
example `() -> only(addprocs(1))` for a replacement on the local machine.
Without `addworker` the pool shrinks when a worker dies.

# Example
The following code will start up 3 processes with Oscar,
run a parallel computation over each element in an array and
then shutdown the processes.
```
results = oscar_worker_pool(3; exeflags="--heap-size-hint=8G") do wp
  Qxy, (x, y) = QQ[:x, :y]
  pmap(z -> z^2, wp, [x^2, x*y, x*y^2])
end
```
"""
oscar_worker_pool(n::Int, init_expr::Expr=:(); kw...) = OscarWorkerPool(n, init_expr; kw...)
oscar_worker_pool(manager::ClusterManager, init_expr::Expr=:(); kw...) = OscarWorkerPool(manager, init_expr; kw...)

function oscar_worker_pool(f::Function, args...; kw...)
  wp = OscarWorkerPool(args...; kw...)
  local results
  try
    results = f(wp)
  finally
    close!(wp)
  end
  yield()
  return results
end

### The following implements the `AbstractWorkerPool` interface;
### `take!`, `put!`, `length` and `isready` are the defaults of Distributed.

# Add a new worker to the pool
function push!(wp::OscarWorkerPool, id::Int)
  # Make sure the node is running Oscar
  remotecall_eval(Main, id, :(using Oscar))
  remotecall_eval(Main, id, wp.init_expr)
  push!(wp.workers, id)
  put!(wp.channel, id)
  return wp
end

### end of implementation interface `AbstractWorkerPool`

### extra functionality

workers(wp::OscarWorkerPool) = collect(wp.workers)

function get_channel(wp::OscarWorkerPool, id::Int; channel_size::Int=1024)
  chnl = get!(wp.oscar_channels, id) do
    RemoteChannel(()->Channel{Any}(channel_size), id)
  end
end

close!(wp::OscarWorkerPool) = rmprocs(workers(wp)...)

# Forget the lost worker `wid` (dead, or alive but to be killed) and start a replacement.
function replace_worker!(wp::OscarWorkerPool, wid::Int)
  wid in procs() && rmprocs(wid; waitfor=0)
  delete!(wp.workers, wid)
  delete!(wp.oscar_channels, wid)
  isnothing(wp.addworker) && return
  try
    push!(wp, wp.addworker())
  catch e
    @warn "could not start a replacement for worker $wid" exception=(e, catch_backtrace())
  end
end

# Hand back the worker `wid` after a call: to the pool if it is alive, replace it if it died.
function release!(wp::OscarWorkerPool, wid::Int)
  wid in procs() ? put!(wp, wid) : replace_worker!(wp, wid)
end

# extend functionality so that `pmap` works with Oscar stuff
function put_type_and_params(wp::OscarWorkerPool, a::Any)
  asyncmap(workers(wp)) do id
    put_type_and_params(get_channel(wp, id), a)
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
# The idea is to use this special type of worker pool to send the `type_and_params`
# of the arguments to the workers up front. Then hopefully, whenever an
# `OscarWorkerPool` is being used, loads of other functionality like `pmap`
# will run out of the box.
=#
function put_args_params(wp::OscarWorkerPool, wid::Int, args, kwargs)
  for a in args
    put_type_and_params(get_channel(wp, wid), a)
  end
  for a in kwargs
    put_type_and_params(get_channel(wp, wid), a)
  end
end

function remotecall(f::Any, wp::OscarWorkerPool, args...; kwargs...)
  wid = take!(wp)
  put_args_params(wp, wid, args, kwargs)

  # Copied from Distributed.jl/src/workerpool.jl.
  # This puts the worker back to the pool once the future is ready
  # and we do not have to worry about this ourselves.
  local fut
  try
    fut = remotecall(f, wid, args...; kwargs...)
  catch
    release!(wp, wid)
    rethrow()
  end
  t = Threads.@spawn Threads.threadpool() try
    wait(fut)
  catch
  finally
    release!(wp, wid)
  end
  errormonitor(t)

  return fut
end

function remotecall_fetch(f::Any, wp::OscarWorkerPool, args...; kwargs...)
  wid = take!(wp)
  # this is inefficient and can clog up channels
  # TODO only send params when necessary
  # and/or figure out how to clear channels
  try
    put_args_params(wp, wid, args, kwargs)
    return remotecall_fetch(f, wid, args...; kwargs...)
  finally
    release!(wp, wid)
  end
end

@doc raw"""
    TimeoutException

Thrown by [`remotecall_with_timeout`](@ref) when the computation did not finish in time.
"""
struct TimeoutException <: Exception
  timeout::Float64
end

Base.showerror(io::IO, e::TimeoutException) = print(io, "computation did not finish within $(e.timeout) seconds")

@doc raw"""
    remotecall_with_timeout(f, wp::OscarWorkerPool, timeout::Real, args...; kwargs...)

Like `remotecall_fetch(f, wp, args...; kwargs...)`, but if the result is not
there after `timeout` seconds, the worker is killed, replaced by a new one (see
[`oscar_worker_pool`](@ref)) and a [`TimeoutException`](@ref) is thrown.
Killing the worker is the only reliable way to stop a computation that is stuck
in a C library.
"""
function remotecall_with_timeout(f::Any, wp::OscarWorkerPool, timeout::Real, args...; kwargs...)
  wid = take!(wp)
  local t
  try
    put_args_params(wp, wid, args, kwargs)
    t = @async remotecall_fetch(f, wid, args...; kwargs...)
  catch
    release!(wp, wid)
    rethrow()
  end
  if timedwait(() -> istaskdone(t), timeout) == :timed_out
    replace_worker!(wp, wid)
    throw(TimeoutException(timeout))
  end
  release!(wp, wid)
  return fetch(t)
end

# this is here so that setting batchsize still works
type_and_params(p::Base.Iterators.Zip) = type_and_params(first(p))

# The pool `pmap` runs on: it shares the free and all workers with `wp`, but its
# calls do not send the `type_and_params` of their arguments again, since `pmap`
# sends those of the whole collection `c` up front; only a worker that joined `wp`
# after that (a replacement for a dead one) gets them on its first call.
# Distributed's `pmap` applies `on_error` on the worker, so it never sees a worker
# that died; with `lost_worker` set (`on_error` without batching) that call returns
# `lost_worker(exception)` instead of ending the whole `pmap`.
struct PmapPool <: AbstractWorkerPool
  channel::Channel{Int}
  workers::Set{Int}
  wp::OscarWorkerPool
  c::Any
  has_params::Set{Int}
  lost_worker::Union{Function, Nothing}
end

function remotecall_fetch(f::Any, p::PmapPool, args...; kwargs...)
  wid = take!(p)
  try
    wid in p.has_params || (put_type_and_params(get_channel(p.wp, wid), p.c); push!(p.has_params, wid))
    return remotecall_fetch(f, wid, args...; kwargs...)
  catch e
    (isnothing(p.lost_worker) || wid in procs()) && rethrow()
    return p.lost_worker(e)
  finally
    release!(p.wp, wid)
  end
end

function pmap(f::Any, wp::OscarWorkerPool, c; on_error=nothing, batch_size=1, kwargs...)
  put_type_and_params(wp, c)
  lost_worker = batch_size == 1 ? on_error : nothing
  p = PmapPool(wp.channel, wp.workers, wp, c, Set(wp.workers), lost_worker)
  return pmap(f, p, c; on_error, batch_size, kwargs...)
end

include("determinants.jl")
export oscar_worker_pool, remotecall_with_timeout, TimeoutException
