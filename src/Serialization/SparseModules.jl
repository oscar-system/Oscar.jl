########################################################################
# Serialization of `FreeMod`s and their elements
########################################################################

@register_serialization_type FreeMod uses_id

type_params(F::T) where {T <: FreeMod} = TypeParams(T, :base_ring=>base_ring(F))

function save_object(s::SerializerState, F::FreeMod)
  save_data_dict(s) do
    #save_object(s, ngens(F), :rank)
    save_object(s, symbols(F), :symbols)
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeMod}, params::Dict)
  R = params[:base_ring]::Ring
  symbs = load_object(s, Vector{Symbol}, :symbols)
  #rk = load_object(s, Int, :rank)
  return FreeMod(R, symbs)
end

@register_serialization_type FreeModElem 

type_params(a::T) where {T<:FreeModElem} = TypeParams(T, parent(a))

save_object(s::SerializerState, v::FreeModElem) = save_object(s, coordinates(v))

function load_object(s::DeserializerState, ::Type{<:FreeModElem}, parent::FreeMod)
  P = base_ring(parent)
  RET = elem_type(P)
  r = load_object(s, SRow, P)
  return FreeModElem(r, parent)
end

########################################################################
# Serialization of `SubquoModule`s and their elements
########################################################################

@register_serialization_type SubquoModule uses_id

type_params(F::T) where {T <: SubquoModule} = TypeParams(T, :base_ring=>base_ring(F))

function save_object(s::SerializerState, M::SubquoModule)
  save_data_dict(s) do
    save_object(s, ambient_free_module(M), :F)
    save_object(s, ambient_representatives_generators(M), :sub)
    save_object(s, relations(M), :quo)
  end
end

function load_object(s::DeserializerState, ::Type{<:SubquoModule}, params::Dict)
  R = params[:base_ring]::Ring
  F = load_object(s, FreeMod{elem_type(R)}, :F)
  I = SubModuleOfFreeModule(F, load_object(s, Vector{FreeModElem{elem_type(R)}}, F, :sub))
  N = SubModuleOfFreeModule(F, load_object(s, Vector{FreeModElem{elem_type(R)}}, F, :quo))
  return SubquoModule(I, N)
end

@register_serialization_type SubquoModuleElem 

type_params(a::T) where {T<:SubquoModuleElem} = TypeParams(T, parent(a))

save_object(s::SerializerState, v::SubquoModuleElem) = save_object(s, coordinates(v))

function load_object(s::DeserializerState, ::Type{<:SubquoModuleElem}, parent::SubquoModule)
  P = base_ring(parent)
  RET = elem_type(P)
  r = load_object(s, SRow, P)
  return SubquoModuleElem(r, parent)
end

