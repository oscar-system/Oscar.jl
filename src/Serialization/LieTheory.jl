using Oscar: _vec, set_root_system_type!

###############################################################################
#
#   Root systems
#
###############################################################################

@register_serialization_type RootSystem uses_id

function save_object(s::SerializerState, R::RootSystem)
  save_data_dict(s) do
    save_object(s, cartan_matrix(R), :cartan_matrix)
    if has_root_system_type(R)
      type, type_ordering = root_system_type_with_ordering(R)
      save_object(s, type, :type)
      if !issorted(type_ordering) # don't save if it's the default
        save_object(s, type_ordering, :type_ordering)
      end
    end
  end
end

function load_object(s::DeserializerState, ::Type{RootSystem})
  cm = load_object(s, Matrix{Int}, :cartan_matrix)
  R = root_system(cm; check=false, detect_type=false)
  if haskey(s, :type)
    type = load_object(s, Vector{Tuple{Symbol,Int}}, :type)
    if haskey(s, :type_ordering)
      type_ordering = load_object(s, Vector{Int}, :type_ordering)
      set_root_system_type!(R, type, type_ordering)
    else
      set_root_system_type!(R, type)
    end
  end
  return R
end

@register_serialization_type RootSpaceElem
@register_serialization_type DualRootSpaceElem

type_params(e::T) where T <: Union{RootSpaceElem, DualRootSpaceElem} = TypeParams(T, root_system(e))

function save_object(s::SerializerState, r::Union{RootSpaceElem,DualRootSpaceElem})
  save_object(s, _vec(coefficients(r)))
end

function load_object(
  s::DeserializerState, T::Type{<:Union{RootSpaceElem,DualRootSpaceElem}}, R::RootSystem
)
  return T(R, load_object(s, Vector{QQFieldElem}))
end

###############################################################################
#
#   Weight lattices
#
###############################################################################

@register_serialization_type WeightLattice uses_id

type_params(P::WeightLattice) = TypeParams(WeightLattice, root_system(P))

function save_object(s::SerializerState, P::WeightLattice)
  save_data_dict(s) do
    # no data required but we leave this function here to generate a valid json
    # and leave root for possible future attrs
  end
end

function load_object(s::DeserializerState, ::Type{WeightLattice}, R::RootSystem)
  return weight_lattice(R)
end

@register_serialization_type WeightLatticeElem

type_params(w::WeightLatticeElem) = TypeParams(WeightLatticeElem, parent(w))

function save_object(s::SerializerState, w::WeightLatticeElem)
  save_object(s, _vec(coefficients(w)))
end

function load_object(s::DeserializerState, ::Type{WeightLatticeElem}, P::WeightLattice)
  return WeightLatticeElem(P, load_object(s, Vector{ZZRingElem}))
end

###############################################################################
#
#   Weyl groups
#
###############################################################################

@register_serialization_type WeylGroup uses_id

type_params(W::WeylGroup) = TypeParams(WeylGroup, root_system(W))

function save_object(s::SerializerState, W::WeylGroup)
  save_data_dict(s) do
    # no data required but we leave this function here to generate a valid json
    # and leave root for possible future attrs
  end
end

function load_object(s::DeserializerState, ::Type{WeylGroup}, R::RootSystem)
  return weyl_group(R)
end

@register_serialization_type WeylGroupElem

type_params(w::WeylGroupElem) = TypeParams(WeylGroupElem, parent(w))

function save_object(s::SerializerState, x::WeylGroupElem)
  save_object(s, word(x))
end

function load_object(s::DeserializerState, ::Type{WeylGroupElem}, W::WeylGroup)
  return W(load_object(s, Vector{UInt8}); normalize=false)
end
