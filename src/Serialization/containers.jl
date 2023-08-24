################################################################################
# Helper functions
get_nested_entry(x::Any) = x
get_nested_entry(v::AbstractArray) = get_nested_entry(v[1])

################################################################################
# Saving and loading vectors

@registerSerializationType(Vector)
type_needs_params(::Type{<:Vector})  = true
const MatVecType{T} = Union{Matrix{T}, Vector{T}}

function save_type_params(s::SerializerState, obj::S) where {T, S <:MatVecType{T}}
  data_dict(s) do
    save_object(s, encode_type(S), :name)
    if type_needs_params(T)
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
  data_array(s) do
    for elem in x
      save_object(s, elem)
    end
  end
end

function load_object(s::DeserializerState, ::Type{<: Vector},
                     v::Vector, params::Type)
  loaded_v = [load_object(s, params, x) for x in v]
  return loaded_v
end

function load_object(s::DeserializerState, ::Type{<: Vector},
                                 v::Vector, params::Tuple)
  T = params[1]
  return [load_object(s, T, x, params[2]) for x in v]
end

function load_object(s::DeserializerState, ::Type{<: Vector},
                     v::Vector, params::Ring)
  T = elem_type(params)
  return [load_object(s, T, x, params) for x in v]
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

function load_object(s::DeserializerState, ::Type{<:Tuple},
                                 v::Vector{Any}, params::Vector)
  return Tuple(
    params[i] isa Type ? load_object(s, params[i], v[i]) :
      load_object(s, params[i][1], v[i], params[i][2]) for i in 1:length(v)
  )
end

################################################################################
# Saving and loading NamedTuple
@registerSerializationType(NamedTuple)
type_needs_params(::Type{<:NamedTuple}) = true

function save_type_params(s::SerializerState, obj::T) where T <: NamedTuple
  data_dict(s) do
    save_object(s, encode_type(NamedTuple), :name)
    s.key = :params
    data_dict(s) do
      s.key = :tuple_params
      data_array(s) do 
        for (i, value) in enumerate(values(obj))
          U = fieldtype(T, i)
          if type_needs_params(U)
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
@registerSerializationType(Matrix)
type_needs_params(::Type{<:Matrix}) = true
  
function save_object(s::SerializerState, mat::Matrix)
  m, n = size(mat)
  data_array(s) do
    for i in 1:m
      save_object(s, [mat[i, j] for j in 1:n])
    end
  end
end

function load_object(s::DeserializerState, ::Type{<:Matrix},
                                 entries::Vector, params::Type)
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
