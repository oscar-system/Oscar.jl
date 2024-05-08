# Add your new types, functions, and methods here.
struct Semigraphoid
  groundset::AbstractVector
  statements::Vector{Tuple{Vector{Int}, Vector{Int}}}
end

function semigraphoid(GS::AbstractVector, S::Vector{Tuple{Vector{Int}, Vector{Int}}})
  return Semigraphoid(GS, S)
end


@register_serialization_type Semigraphoid

function save_object(s::SerializerState, sm_gp::Semigraphoid)
  save_data_dict(s) do
    save_object(s, sm_gp.groundset, :groundset)
    save_object(s, sm_gp.statements, :statements)
  end
end

function load_object(s::DeserializerState, ::Type{Semigraphoid})
  GS = load_object(s, Vector{Int}, :groundset)
  S = load_object(s, Vector, (Tuple, [(Vector, Int), (Vector, Int)]), :statements)
  return semigraphoid(GS, S)
end






























