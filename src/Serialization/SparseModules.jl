#=
# TODO: Serialization of sparse rows makes trouble with the code below. I don't know why.
@register_serialization_type SRow 

type_params(v::T) where {T<:SRow} = TypeParams(T, base_ring(v))

function save_object(s::SerializerState, v::SRow)
  save_object(s, v.pos)
  save_object(s, v.values)
end

function load_object(s::DeserializerState, ::Type{<:SRow}, parent_ring::Ring)
  pos = load_object(s, Vector{Int}, 1)
  vals = load_object(s, Vector{elem_type(parent_ring)}, parent_ring, 2)
  return sparse_row(parent_ring, pos, vals)
end
=#

########################################################################
# Serialization of `FreeMod`s and their elements
########################################################################

@register_serialization_type FreeMod uses_id

type_params(F::T) where {T <: FreeMod} = TypeParams(T, base_ring(F))

function save_object(s::SerializerState, F::FreeMod)
  save_data_dict(s) do
    #save_object(s, ngens(F), :rank)
    save_object(s, symbols(F), :symbols)
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeMod}, params::Ring)
  R = params[:base_ring]::Ring
  symbs = load_object(s, Vector{Symbol}, :symbols)
  #rk = load_object(s, Int, :rank)
  return FreeMod(R, symbs)
end

@register_serialization_type FreeModElem 

type_params(a::T) where {T<:FreeModElem} = TypeParams(T, parent(a))

function save_object(s::SerializerState, v::FreeModElem)
  r = coordinates(v)
  save_data_dict(s) do
    save_object(s, r.pos, :positions)
    save_object(s, r.values, :values)
  end
end

function load_object(s::DeserializerState, ::Type{<:FreeModElem}, parent::FreeMod)
  P = base_ring(parent)
  RET = elem_type(P)
  pos = load_object(s, Vector{Int}, :positions)
  vals = load_object(s, Vector{RET}, P, :values)
  r = sparse_row(P, pos, vals)
  return FreeModElem(r, parent)
end

########################################################################
# Serialization of `SubquoModule`s and their elements
########################################################################

@register_serialization_type SubquoModule uses_id

type_params(M::T) where {T <: SubquoModule} = TypeParams(T, ambient_free_module(M))

function save_object(s::SerializerState, M::SubquoModule)
  save_data_dict(s) do
    save_object(s, ambient_representatives_generators(M), :sub)
    save_object(s, relations(M), :quo)
  end
end

function load_object(s::DeserializerState, ::Type{<:SubquoModule}, params::FreeMod)
  F = params
  R = base_ring(F)
  I = SubModuleOfFreeModule(F, load_object(s, Vector{FreeModElem{elem_type(R)}}, F, :sub))
  N = SubModuleOfFreeModule(F, load_object(s, Vector{FreeModElem{elem_type(R)}}, F, :quo))
  return SubquoModule(I, N)
end

@register_serialization_type SubquoModuleElem 

type_params(a::T) where {T<:SubquoModuleElem} = TypeParams(T, parent(a))

function save_object(s::SerializerState, v::SubquoModuleElem)
  r = coordinates(v)
  save_data_dict(s) do
    save_object(s, r.pos, :positions)
    save_object(s, r.values, :values)
  end
end

function load_object(s::DeserializerState, ::Type{<:SubquoModuleElem}, parent::SubquoModule)
  P = base_ring(parent)
  RET = elem_type(P)
  pos = load_object(s, Vector{Int}, :positions)
  vals = load_object(s, Vector{RET}, P, :values)
  r = sparse_row(P, pos, vals)
  return SubquoModuleElem(r, parent)
end

