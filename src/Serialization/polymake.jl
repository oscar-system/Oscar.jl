# Loading polymake objects into the appropriate types.


##############################################################################
"""
    load_from_polymake(::Type{Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Directed, Undirected}}

Load a graph stored in JSON format, given the filename as input.
"""
function load_from_polymake(::Type{Graph{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Directed, Undirected}}
  polymake_object = Polymake.call_function(:common, :deserialize_json_string, JSON3.write(jsondict))
  return Graph{T}(polymake_object)
end

function load_from_polymake(::Type{PhylogeneticTree{T}}, jsondict::Dict{Symbol, Any}) where {T <: Union{Float64, Int, QQFieldElem}}
  polymake_object = Polymake.call_function(:common, :deserialize_json_string, JSON3.write(jsondict))
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
  inner_object = Polymake.call_function(:common, :deserialize_json_string, JSON3.write(jsondict))
  if T <: PolyhedralObject{Float64}
    return T(inner_object, AbstractAlgebra.Floats{Float64}())
  end
  try
    return T(inner_object)
  catch e
    error("Unsupported object type $T for loading polymake object")
  end
end

# Distinguish between the various polymake datatypes.
function load_from_polymake(jsondict::Dict{Symbol, Any})
  typename = jsondict[:_type]
  if haskey(polymake2OscarTypes, typename)
    oscar_type = polymake2OscarTypes[typename]
    return load_from_polymake(oscar_type, jsondict)
  else 
    # We just try to default to something from Polymake.jl
    deserialized = Polymake.call_function(:common, :deserialize_json_string, JSON3.write(jsondict))
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

#_pmdata_for_oscar(bo::Polymake.BigObject, coeff::Field) = _bigobject_to_dict(bo, coeff)

_pmdata_for_oscar(::Nothing, coeff::Field) = nothing
_pmdata_for_oscar(v::Union{Bool,Int64,Float64,String}, coeff::Field) = v
if Polymake.CxxWrap.CxxLong != Int64
  _pmdata_for_oscar(i::Polymake.CxxWrap.CxxLong, coeff::Field) = Int64(i)
end

_pmdata_for_oscar(im::IncidenceMatrix, coeff::Field) = im

_pmdata_for_oscar(g::Polymake.Graph{T}, coeff::Field) where T = Graph{T}(g)

_pmdata_for_oscar(m::Polymake.Matrix, coeff::Field) = matrix(coeff, m)
_pmdata_for_oscar(m::Polymake.Matrix{Polymake.to_cxx_type(Int)}, coeff::Field) = Matrix{Int}(m)
_pmdata_for_oscar(m::Polymake.Matrix{<:Polymake.Integer}, coeff::Field) = matrix(ZZ, m)
_pmdata_for_oscar(m::Polymake.Matrix{<:Polymake.Rational}, coeff::Field) = matrix(QQ, m)

_pmdata_for_oscar(m::Polymake.SparseMatrix, coeff::Field) = _pmdata_for_oscar(Polymake.common.dense(m), coeff)

_pmdata_for_oscar(v::Polymake.Vector, coeff::Field) = collect(elem_type(coeff), map(coeff, v))
_pmdata_for_oscar(v::Polymake.Vector{Polymake.to_cxx_type(Int)}, coeff::Field) = Vector{Int}(v)
_pmdata_for_oscar(v::Polymake.Vector{<:Polymake.Integer}, coeff::Field) = collect(ZZRingElem, map(ZZ, v))
_pmdata_for_oscar(v::Polymake.Vector{<:Polymake.Rational}, coeff::Field) = collect(QQFieldElem, map(QQ, v))

_pmdata_for_oscar(v::Polymake.SparseVector, coeff::Field) = _pmdata_for_oscar(Polymake.common.dense(v), coeff)

_pmdata_for_oscar(nm::Polymake.NodeMap, coeff::Field) = _pmdata_for_oscar(Polymake.Array(nm), coeff)
_pmdata_for_oscar(bd::Polymake.BasicDecoration, coeff::Field) = (_pmdata_for_oscar(Polymake.decoration_face(bd), coeff), _pmdata_for_oscar(Polymake.decoration_rank(bd), coeff))

_pmdata_for_oscar(s::Polymake.Integer, coeff::Field) = ZZ(s)
_pmdata_for_oscar(s::Polymake.Rational, coeff::Field) = QQ(s)
_pmdata_for_oscar(s::Polymake.OscarNumber, coeff::Field) = coeff(s)

_convert_pm_minormax(::Type{Polymake.Max}) = typeof(max)
_convert_pm_minormax(::Type{Polymake.Min}) = typeof(min)
_pmdata_for_oscar(s::Polymake.TropicalNumber{A}, coeff::Field) where A = tropical_semiring(_convert_pm_minormax(A))(s)

_pmdata_for_oscar(s::Polymake.CxxWrap.StdString, coeff::Field) = String(s)

_pmdata_for_oscar(a::Polymake.Array, coeff::Field) = [_pmdata_for_oscar(e, coeff) for e in a]
_pmdata_for_oscar(a::Polymake.Array{T}, coeff::Field) where T <: Union{Polymake.Matrix, Polymake.Vector} = Tuple(_pmdata_for_oscar.(a, Ref(coeff)))

_pmdata_for_oscar(s::Polymake.Set, coeff::Field) = Set(_pmdata_for_oscar(e, coeff) for e in s)


function _bigobject_to_dict(bo::Polymake.BigObject, coeff::Field, parent_key::String="")
  data = Dict{Symbol,Any}()
  for pname in Polymake.list_properties(bo)
    p = Polymake.give(bo, pname)
    key_str = parent_key == "" ? pname : parent_key * "." * pname
    if p isa Polymake.PropertyValue
      @debug "missing c++ mapping: skipping $pname of type $(Polymake.typeinfo_string(p, true))"
    elseif p isa Polymake.BigObject
       obj = _bigobject_to_dict(p, coeff, key_str)
       merge!(data, obj)
    else
      try
        obj = _pmdata_for_oscar(p, coeff)
        data[Symbol(key_str)] = obj
      catch e
        if e isa MethodError
          @debug "failed to convert $pname of type $(typeof(p)) to Oscar, skipping"
        else
          rethrow(e)
        end
      end
    end
  end
  description = Polymake.getdescription(bo)

  if !isempty(description)
    data[:_description] = String(description)
  end
  data
end

function _polyhedral_object_as_dict(x::Oscar.PolyhedralObjectUnion)
  bo = Oscar.pm_object(x)
  data = _bigobject_to_dict(bo, coefficient_field(x))
  data[:_polymake_type] = Polymake.bigobject_qualifiedname(bo)
  data[:_coeff] = coefficient_field(x)
  return data
end

function _load_bigobject_from_dict!(obj::Polymake.BigObject, dict::Dict, parent_key::String="")
  delay_loading = Tuple{Symbol,Any}[]
  for (k, v) in dict
    # keys of dict are symbols
    k = string(k)
    key_str = parent_key == "" ? k : parent_key * "." * k
    first(k) == '_' && continue

    if v isa Dict
      _load_bigobject_from_dict!(obj, v, key_str)
    else
      pmv = convert(Polymake.PolymakeType, v)
      bot = Polymake.bigobject_type(obj)
      # NodeMaps need extra treatment since the constructor doesn't accept polymake c++ arrays
      # and we can't create a nodemap from scratch without the graph
      # so we convert it to a pure perl array and delay loading until the end of this level
      if pmv isa Polymake.Array && Polymake.bigobject_prop_type(bot, key_str) in ["NodeMap", "EdgeMap"]
        pmv = Polymake.as_perl_array_of_array(pmv)
        push!(delay_loading, (Symbol(key_str), pmv))
      else
        Polymake.take(obj, key_str, pmv)
      end
    end
  end
  for (k, v) in delay_loading
    Polymake.take(obj, string(k), v)
  end
  if haskey(dict, :_description)
    Polymake.setdescription!(obj, dict[:_description])
  end
  return obj
end

function _dict_to_bigobject(dict::Dict{Symbol, Any})
  obj = Polymake.BigObject(Polymake.BigObjectType(dict[:_polymake_type]))
  _load_bigobject_from_dict!(obj, dict)
  return obj
end
