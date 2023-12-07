using Distributed

function set_channels(input::Type{S}, output::Type{T}, params::Type{U};
                      channel_size::Tuple{Int, Int, Int} = (32, 32, 32)) where {S, T, U}
  return (
    RemoteChannel(()->Channel{S}(channel_size[1])),
    [RemoteChannel(() -> Channel{T}(channel_size[2]), w) for w in workers()],
    [RemoteChannel(() -> Channel{U}(channel_size[3]), w) for w in workers()]
  )
end

function put_params(params_channels::Vector{<: RemoteChannel}, params::Any)
  for (i, w) in enumerate(workers())
    put!(params_channels[i], params)
  end

  function take_params(chnnl)
    params = take!(chnnl)
  end

  for (i, w) in enumerate(workers())
    remotecall(take_params, w, params_channels[i]) 
  end
end

  

