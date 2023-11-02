################################################################################
# Helper functions
get_nested_entry(x::Any) = x
get_nested_entry(v::AbstractArray) = get_nested_entry(v[1])

################################################################################
# Saving and loading vectors

@register_serialization_type Vector uses_params

const MatVecType{T} = Union{Matrix{T}, Vector{T}}

function save_type_params(s::SerializerState, obj::S) where {T, S <:MatVecType{T}}
  save_data_dict(s) do
    save_object(s, encode_type(S), :name)
    if serialize_with_params(T)
      if hasmethod(parent, (T,))
        parents = map(parent, obj)
        parents_all_equal = all(map(x -> isequal(first(parents), x), parents))
        @req parents_all_equal "Not all parents of Vector or Matrix entries are the same, consider using a Tuple"
      end
      save_type_params(s, obj[1], :params)
    else
      save_object(s, encode_type(T), :params)
    end
  end
end

function load_type_params(s::DeserializerState, ::Type{<:MatVecType}, str::String)
  return decode_type(str)
end

function load_type_params(s::DeserializerState, ::Type{<:MatVecType}, dict::Dict)
  T = decode_type(dict[:name])
  params = load_type_params(s, T, dict[:params])
  return (T, params)
end

function load_type_params(s::DeserializerState, ::Type{<:MatVecType}, override_params::Any)
  return (elem_type(override_params), override_params)
end

function save_object(s::SerializerState, x::Vector)
  save_data_array(s) do
    for elem in x
      if serialize_with_id(typeof(elem))
        ref = save_as_ref(s, elem)
        save_object(s, ref)
      else
        save_object(s, elem)
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: Vector},
                     v::Vector, params::Type)
  if serialize_with_id(params)
    loaded_v = params[load_ref(s, x) for x in v]
  else
    loaded_v = params[load_object(s, params, x) for x in v]
  end
  return loaded_v
end

# handles nested Vectors
function load_object(s::DeserializerState, ::Type{<: Vector},
                                 v::Vector, params::Tuple)
  T = params[1]
  if isempty(v)
    return T[]
  else
    return [load_object(s, T, x, params[2]) for x in v]
  end
end

function load_object(s::DeserializerState, ::Type{<: Vector},
                     v::Vector, params::Ring)
  T = elem_type(params)
  return T[load_object(s, T, x, params) for x in v]
end

################################################################################
# Saving and loading Tuple
@register_serialization_type Tuple uses_params

function save_type_params(s::SerializerState, tup::T) where T <: Tuple
  save_data_dict(s) do
    save_object(s, encode_type(Tuple), :name)
    n = fieldcount(T)
    save_data_array(s, :params) do 
      for i in 1:n
        U = fieldtype(T, i)
        if serialize_with_params(U)
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
  save_data_array(s) do 
    for entry in obj
      if serialize_with_id(typeof(entry))
        ref = save_as_ref(s, entry)
        save_object(s, ref)
      else
        save_object(s, entry)
      end
    end
  end
end

function get_nested_type(params::Vector)
  if params[2] isa Type
    return params[1]{params[2]}
  end
  nested_type = get_nested_type(params[2])
  return params[1]{nested_type}
end

function get_tuple_type(params::Vector)
  type_vector = Type[]

  for t in params
    if t isa Type
      push!(type_vector, t)
    else
      push!(type_vector, get_nested_type(t))
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:Tuple},
                     v::Vector{Any}, params::Vector)
  entries = []
  for i in 1:length(v)
    if params[i] isa Type
      if serialize_with_id(params[i])
        loaded_obj = load_ref(s, v[i])
      else
        loaded_obj = load_object(s, params[i], v[i])
      end
    else
      loaded_obj = load_object(s, params[i][1], v[i], params[i][2])
    end
    push!(entries, loaded_obj)
  end
  return Tuple(entries)
end

################################################################################
# Saving and loading NamedTuple
@register_serialization_type NamedTuple uses_params

function save_type_params(s::SerializerState, obj::T) where T <: NamedTuple
  save_data_dict(s) do
    save_object(s, encode_type(NamedTuple), :name)
    save_data_dict(s, :params) do
      save_data_array(s, :tuple_params) do 
        for (i, value) in enumerate(values(obj))
          U = fieldtype(T, i)
          if serialize_with_params(U)
            save_type_params(s, value)
          else
            save_object(s, encode_type(U))
          end
        end
      end
      save_object(s, keys(obj), :names)
    end
  end
end

function save_object(s::SerializerState, obj::NamedTuple)
  save_object(s, values(obj))
end

function load_type_params(s::DeserializerState, ::Type{<:NamedTuple}, params::Dict)
  loaded_params = Any[]
  for param in params[:tuple_params]
    if param isa String
      push!(loaded_params, decode_type(param))
    else
      T = decode_type(param[:name])
      push!(loaded_params, (T, load_type_params(s, T, param[:params])))
    end
  end
  return (params[:names], loaded_params)
end

function load_object(s::DeserializerState, ::Type{<: NamedTuple},
                                 v::Vector, params::Tuple)
  keys, tuple_params = params
  tuple = load_object(s, Tuple, v, tuple_params)
  keys = map(Symbol, keys)
  return NamedTuple{Tuple(keys), typeof(tuple)}(tuple)
end

################################################################################
# Saving and loading matrices
@register_serialization_type Matrix uses_params
  
function save_object(s::SerializerState, mat::Matrix)
  m, n = size(mat)
  save_data_array(s) do
    for i in 1:m
      save_object(s, [mat[i, j] for j in 1:n])
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:Matrix},
                     entries::Vector, params::Type)
  if isempty(entries)
    return zero_matrix(parent_type(params)(), 0, 0)
  end
  m = reduce(vcat, [
    permutedims(load_object(s, Vector, v, params)) for v in entries
      ])

  return Matrix{params}(m)
end

function load_object(s::DeserializerState, ::Type{<:Matrix},
                     entries::Vector, params::Tuple)
  m = reduce(vcat, [
    permutedims(load_object(s, Vector, v, params)) for v in entries
      ])
  return Matrix{params[1]}(m)
end
