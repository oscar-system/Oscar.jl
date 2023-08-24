@registerSerializationType(ToricCoveredScheme{QQField})

function save_object(s::SerializerState, tcs::ToricCoveredScheme{QQField})
  save_typed_object(s, tcs.ntv)
end

function load_object(s::DeserializerState, ::Type{ToricCoveredScheme{QQField}}, dict::Dict)
  return ToricCoveredScheme{QQField}(load_typed_object(s, dict))
end

@registerSerializationType(GlobalTateModel)

function save_object(s::SerializerState, gtm::GlobalTateModel)
  data_dict(s) do
    s.key = :tate_as
    data_array(s) do
      # I am not sure what the type of these entries are? here you might not want to
      # serialize as a Vector but use another data_array(s) do ... 
      save_object(s, [gtm.tate_a1, gtm.tate_a2, gtm.tate_a3, gtm.tate_a4, gtm.tate_a6])
      save_object(s, gtm.tate_polynomial, :tate_polynomial)

      # this property can probably be moved to type parameters
      save_object(s, gtm.base_space, :base_space) 
    end
end
