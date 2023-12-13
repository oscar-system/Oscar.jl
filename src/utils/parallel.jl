using Distributed

"""
  params_channels(params::Type{U}, proc_ids::Vector{Int} = workers(), channel_size::Int = 32) where U

Sets the parameter channel to be used for sending intermediate parameters used during serialization.
Since the IPCserializer doesn't serialize it's references, it is up to the user to guarantee that
the objects corresponding to these references exist on each process. 
"""
function params_channels(params::Type{U}, proc_ids::Vector{Int} = workers(),
                         channel_size::Int = 32) where U
  return [RemoteChannel(() -> Channel{U}(channel_size), p) for p in proc_ids]
end


"""
  put_params(params_channels::Vector{<: RemoteChannel}, params::Any, proc_ids::Vector{Int} = workers())

Uses the parameter channel and sends `params` to each process with id in `proc_ids`
"""
function put_params(params_channels::Vector{<: RemoteChannel}, params::Any,
                    proc_ids::Vector{Int} = workers())
  for (i, w) in enumerate(proc_ids)
    put!(params_channels[i], params)
  end

  function take_params(chnnl)
    params = take!(chnnl)
  end

  for (i, w) in enumerate(proc_ids)
    remotecall(take_params, w, params_channels[i]) 
  end
end

  

