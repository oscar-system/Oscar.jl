# Loading polymake objects into the appropriate types.


##############################################################################
"""
    load_from_polymake(::Type{Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Directed, Undirected}}

Load a graph stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Directed, Undirected}}
  polymake_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
  return Graph{T}(polymake_object)
end

function load_from_polymake(::Type{PhylogeneticTree{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Float64, Int, QQFieldElem}}
  polymake_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
  return PhylogeneticTree{T}(polymake_object)
end

const polymake2OscarTypes = Dict{String, Type}([
  "polytope::Cone<Rational>" => Cone{QQFieldElem},
  "polytope::Polytope<Rational>" => Polyhedron{QQFieldElem},
  "fan::PolyhedralFan<Rational>" => PolyhedralFan{QQFieldElem},
  "fan::PolyhedralComplex<Rational>" => PolyhedralComplex{QQFieldElem},
  "fan::SubdivisionOfPoints<Rational>" => SubdivisionOfPoints{QQFieldElem},
  "topaz::SimplicialComplex" => SimplicialComplex,
  "common::GraphAdjacency<Undirected>" => Graph{Undirected},
  "common::GraphAdjacency<Directed>" => Graph{Directed},
  "graph::PhylogeneticTree<Rational>" => PhylogeneticTree{QQFieldElem},
  "graph::PhylogeneticTree<Float>" => PhylogeneticTree{Float64}
])

@register_serialization_type Polymake.BigObjectAllocated "Polymake.BigObject" uses_id

function load_from_polymake(::Type{T}, jsondict::Dict{Symbol, Any}) where {
  T<:Union{Cone{<:scalar_types}, Polyhedron{<:scalar_types}, PolyhedralFan{<:scalar_types}, 
           PolyhedralComplex{<:scalar_types}, SubdivisionOfPoints{<:scalar_types}, SimplicialComplex}}
  
  inner_object = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
  return T(inner_object)
end


# Distinguish between the various polymake datatypes.
function load_from_polymake(jsondict::Dict{Symbol, Any})
  typename = jsondict[:_type]
  if haskey(polymake2OscarTypes, typename)
    oscar_type = polymake2OscarTypes[typename]
    return load_from_polymake(oscar_type, jsondict)
  else 
    # We just try to default to something from Polymake.jl
    deserialized = Polymake.call_function(:common, :deserialize_json_string, json(jsondict))
    if !isa(deserialized, Polymake.BigObject)
      @warn "No function for converting the deserialized Polymake type to Oscar type: $(typeof(deserialized))"
      return deserialized
    end
    
    if Polymake.type_name(deserialized) == "Ideal"
      return convert(MPolyIdeal{QQMPolyRingElem}, deserialized)
    else
      @warn "No function for converting the deserialized Polymake type to Oscar type: $(typeof(deserialized))"
      return deserialized
    end
  end
end

