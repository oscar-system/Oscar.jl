################################################################################
# Helper functions
get_nested_entry(x::Any) = x
get_nested_entry(v::AbstractArray) = get_nested_entry(v[1])

################################################################################
# Saving and loading vectors

@registerSerializationType(Vector)
type_needs_params(::Type{<:Vector})  = true

function save_type_params(s::SerializerState, obj::Vector{T}) where T
  data_dict(s) do
    save_object(s, encode_type(Vector), :name)
    if type_needs_params(T)
        save_type_params(s, obj[1], :params)
    else
      save_object(s, encode_type(T), :params)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:Vector}, str::String)
  return decode_type(str)
end

function load_type_params(s::DeserializerState, ::Type{<:Vector}, dict::Dict)
  T = decode_type(dict[:name])
  params = load_type_params(s, T, dict[:params])
  return (T, params)
end

function load_type_params(s::DeserializerState, ::Type{<:Vector}, override_params::Any)
  return (elem_type(override_params), override_params)
end

function load_object_with_params(s::DeserializerState, ::Type{<: Vector},
                                 v::Vector{Any}, params::Type)
  loaded_v = params[load_object(s, params, x) for x in v]
  return loaded_v
end

function load_object_with_params(s::DeserializerState, ::Type{<: Vector},
                                 v::Vector{Any}, params::Tuple)
  T = params[1]
  return [load_object_with_params(s, T, x, params[2]) for x in v]
end

################################################################################
# Saving and loading Tuple
@registerSerializationType(Tuple)
type_needs_params(::Type{<:Tuple}) = true

function save_type_params(s::SerializerState, tup::T) where T <: Tuple
  data_dict(s) do
    save_object(s, encode_type(Tuple), :name)
    n = fieldcount(T)
    s.key = :params
    data_array(s) do 
      for i in 1:n
        U = fieldtype(T, i)
        if type_needs_params(U)
          save_type_params(s, tup[i])
        else
          save_object(s, encode_type(U))
        end
      end
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:Tuple}, params::Vector)
  loaded_params = Any[]
  for param in params
    if param isa String
      push!(loaded_params, decode_type(param))
    else
      T = decode_type(param[:name])
      push!(loaded_params, (T, load_type_params(s, T, param[:params])))
    end
  end
  return loaded_params
end

function save_object(s::SerializerState, obj::Tuple)
  data_array(s) do 
    for entry in obj
      save_object(s, entry)
    end
  end
end

function load_object_with_params(s::DeserializerState, ::Type{<:Tuple},
                                 v::Vector{Any}, params::Vector)
  return Tuple(
    params[i] isa Type ? load_object(s, params[i], v[i]) :
      load_object_with_params(s, params[i][1], v[i], params[i][2]) for i in 1:length(v)
  )
end


function load_internal(s::DeserializerState, T::Type{<:Tuple}, dict::Dict)
  field_types = map(decode_type, dict[:field_types])
  n = length(field_types)
  content = dict[:content]
  @assert length(content) == n  "Wrong length of tuple, data may be corrupted."
  return T(load_type_dispatch(s, field_types[i], content[i]) for i in 1:n)
end

################################################################################
# Saving and loading NamedTuple
@registerSerializationType(NamedTuple)

function save_internal(s::SerializerState, n_tup:: NamedTuple)
  return Dict(
    :keys => save_type_dispatch(s, keys(n_tup)),
    :content => save_type_dispatch(s, values(n_tup))
  )
end

function load_internal(s::DeserializerState, ::Type{<:NamedTuple}, dict::Dict)
  tup = load_unknown_type(s, dict[:content])
  keys = load_type_dispatch(s, Tuple, dict[:keys])
  return NamedTuple{keys}(tup)
end


################################################################################
# Saving and loading matrices

@registerSerializationType(Matrix)

function save_internal(s::SerializerState, mat::Matrix{T}) where T
  m, n = size(mat)
  return Dict(
    :matrix => [save_type_dispatch(s, [mat[i, j] for j in 1:n]) for i in 1:m]
  )
end

# deserialize with specific content type
function load_internal(s::DeserializerState, ::Type{Matrix{T}}, dict::Dict) where T
  x = dict[:matrix]
  y = reduce(vcat, [permutedims(load_type_dispatch(s, Vector{T}, x[i])) for i in 1:length(x)])
  return Matrix{T}(y)
end

# deserialize without specific content type
function load_internal(s::DeserializerState, ::Type{Matrix}, dict::Dict)
  x = dict[:matrix]
  y = reduce(vcat, [permutedims(load_type_dispatch(s, Vector, x[i])) for i in 1:length(x)])
  return Matrix(y)
end

# deserialize without specific content type
function load_internal_with_parent(s::DeserializerState,
                                   ::Type{Matrix}, dict::Dict, parent)
  x = dict[:matrix]
  y = reduce(vcat, [
    permutedims(load_type_dispatch(s, Vector, x[i], parent=parent)) for i in 1:length(x)
      ])
  return Matrix(y)
end
