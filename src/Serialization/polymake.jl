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

_pmdata_for_oscar(bo::Polymake.BigObject, coeff::Union{Field,Nothing}) = _bigobject_to_dict(bo, coeff)

_pmdata_for_oscar(::Nothing, coeff::Union{Field,Nothing}) = nothing
_pmdata_for_oscar(v::Union{Bool,Int64,Float64,String}, coeff::Union{Field,Nothing}) = v
if Polymake.CxxWrap.CxxLong != Int64
  _pmdata_for_oscar(i::Polymake.CxxWrap.CxxLong, coeff::Union{Field,Nothing}) = Int64(i)
end

_pmdata_for_oscar(im::IncidenceMatrix, coeff::Union{Field,Nothing}) = im

_pmdata_for_oscar(g::Polymake.Graph{T}, coeff::Union{Field,Nothing}) where T = Graph{T}(g)

_pmdata_for_oscar(m::Polymake.Matrix, coeff::Union{Field,Nothing}) = matrix(coeff, m)
_pmdata_for_oscar(m::Polymake.Matrix{<:Polymake.Integer}, coeff::Union{Field,Nothing}) = matrix(ZZ, m)
_pmdata_for_oscar(m::Polymake.Matrix{<:Polymake.Rational}, coeff::Union{Field,Nothing}) = matrix(QQ, m)

_pmdata_for_oscar(m::Polymake.SparseMatrix, coeff::Union{Field,Nothing}) = _pmdata_for_oscar(Polymake.common.dense(m), coeff)

_pmdata_for_oscar(v::Polymake.Vector, coeff::Union{Field,Nothing}) = collect(tt, map(F, v))
_pmdata_for_oscar(v::Polymake.Vector{<:Polymake.Integer}, coeff::Union{Field,Nothing}) = collect(ZZRingElem, map(ZZ, v))
_pmdata_for_oscar(v::Polymake.Vector{<:Polymake.Rational}, coeff::Union{Field,Nothing}) = collect(QQFieldElem, map(QQ, v))

_pmdata_for_oscar(v::Polymake.SparseVector, coeff::Union{Field,Nothing}) = _pmdata_for_oscar(Polymake.common.dense(v), coeff)

_pmdata_for_oscar(s::Polymake.Integer, coeff::Union{Field,Nothing}) = ZZ(s)
_pmdata_for_oscar(s::Polymake.Rational, coeff::Union{Field,Nothing}) = QQ(s)
_pmdata_for_oscar(s::Polymake.OscarNumber, coeff::Union{Field,Nothing}) = coeff(s)

_pmdata_for_oscar(s::Polymake.CxxWrap.StdString, coeff::Union{Field,Nothing}) = String(s)

_pmdata_for_oscar(a::Polymake.Array, coeff::Union{Field,Nothing}) = [_pmdata_for_oscar(e, coeff) for e in a]
_pmdata_for_oscar(s::Polymake.Set, coeff::Union{Field,Nothing}) = Set(_pmdata_for_oscar(e, coeff) for e in s)


function _bigobject_to_dict(bo::Polymake.BigObject, coeff::Union{Field,Nothing})
  data = Dict{String,Any}()
  data["_type"] = Polymake.bigobject_qualifiedname(bo)
  data["_coeff"] = coeff
  for pname in Polymake.list_properties(bo)
    p = Polymake.give(bo, pname)
    if p isa Polymake.PropertyValue
      @warn "missing c++ mapping: skipping $pname of type $(Polymake.typeinfo_string(p, true))"
    else
      try
        obj = _pmdata_for_oscar(p, coeff)
        data[pname] = obj
      catch MethodError
        @warn "failed to convert $pname of type $(typeof(p)) to Oscar, skipping"
      end
    end
  end
  data
end

function _polyhedral_object_as_dict(x::Oscar.PolyhedralObjectUnion)
  return _bigobject_to_dict(Oscar.pm_object(x), coefficient_field(x))
end

function _load_bigobject_from_dict!(obj::Polymake.BigObject, dict::Dict, parent_key::String="")
  for (k, v) in dict
    key_str = parent_key == "" ? k : parent_key * "." * k

    k == "_type" && continue
    k == "_coeff" && continue

    if v isa Dict
      _load_bigobject_from_dict!(obj, v, key_str)
    else
      Polymake.take(obj, key_str, convert(Polymake.PolymakeType, v))
    end
  end
end

function _dict_to_bigobject(dict::Dict{String, Any})
  obj = Polymake.BigObject(Polymake.BigObjectType(dict["_type"]))
  _load_bigobject_from_dict!(obj, dict)
  return obj
end
