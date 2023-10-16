## TODO
# it might be beneficial to eventually implement a save_type_params / load_type_params
# especially if one wants to serialize vectors of these objects
# also ZZLatWithIsom currently always serialized with ambient space, this may not alwyas be
# desired

@register_serialization_type QuadSpaceWithIsom

function save_object(s::SerializerState, QS::QuadSpaceWithIsom)
  save_data_dict(s) do
    save_typed_object(s, space(QS), :quad_space)
    save_typed_object(s, isometry(QS), :isom)

    save_object(s, order_of_isometry(QS), :order)
  end
end

function load_object(s::DeserializerState, ::Type{QuadSpaceWithIsom}, dict::Dict)
  quad_space = load_typed_object(s, dict[:quad_space])
  isom = load_typed_object(s, dict[:isom])

  # not quite sure how to deal with IntExt/PosInf yet..
  # we could add it to the basic type section
  # then we could use
  # load_object(s, IntExt, dict[:order]) 
  n = load_object(s, Int, dict[:order])
  
  return QuadSpaceWithIsom(quad_space, isom, n)
end

@register_serialization_type ZZLatWithIsom

function save_object(s::SerializerState, L::ZZLatWithIsom)
  save_data_dict(s) do
    save_typed_object(s, ambient_space(L), :ambient_space)
    save_typed_object(s, lattice(L), :lattice)
  end
end

function load_object(s::DeserializerState, ::Type{ZZLatWithIsom}, dict::Dict)
  quad_space = load_typed_object(s, dict[:ambient_space])
  lat = load_typed_object(s, dict[:lattice])

  return ZZLatWithIsom(quad_space, lat, isometry(quad_space), order_of_isometry(quad_space))
end
