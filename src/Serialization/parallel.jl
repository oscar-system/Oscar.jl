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

function put_type_params(channel::RemoteChannel, tps::Tuple{Vararg{TypeParams}})
  for tp in tps
    put_type_params(channel, params(tp))
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