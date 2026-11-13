################################################################################
# The mrdi namespace
#
# Types that are agreed upon between the systems implementing the mrdi format
# form the `mrdi` namespace. It is the default namespace of a file whose header
# lists `mrdi`, so these types are written without a namespace tag. All other
# types are tagged with `"_ns": "Oscar"`.

const mrdi_namespace_url = "https://www.mardi4nfdi.de/mrdi"
const mrdi_namespace_version = v"1.0.0"

const mrdi_types = Set{String}([
  "Boolean",
  "String",
  "Symbol",
  "Integers",
  "Integer",
  "Rationals",
  "Rational",
  "HomogeneousList",
  "InhomogeneousList",
  "Set",
  "FiniteField",
  "FiniteFieldElement",
  "PolynomialRing",
  "Polynomial",
])

is_mrdi_type(T::Type) = encode_type(T) in mrdi_types

# writes the name of `T` and, if needed, its namespace into the current dict
function save_type_name(s::SerializerState, T::Type)
  save_object(s, encode_type(T), :name)
  is_mrdi_type(T) || save_object(s, "Oscar", :_ns)
end

# writes `T` as a type without params, at `key` or at the key already set
function save_type_name(s::SerializerState, T::Type, key::Union{Symbol, Nothing})
  isnothing(key) || set_key(s, key)
  if is_mrdi_type(T)
    save_object(s, encode_type(T))
  else
    save_data_dict(s) do
      save_type_name(s, T)
    end
  end
end

# whether the current node is a type without params, that is a type name
# or a dict holding only the name and the namespace of a type
function is_type_without_params(s::DeserializerState)
  node_is_string(s) && return true
  return node_is_object(s) && haskey(s, :name) && haskey(s, :_ns) &&
    !haskey(s, :params) && !haskey(s, :_instance)
end
