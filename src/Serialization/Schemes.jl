### affine schemes

@register_serialization_type AffineScheme uses_id

function save_object(s::SerializerState, X::AffineScheme)
  save_data_dict(s) do
    save_typed_object(s, OO(X), :coordinate_ring)
  end
end

function load_object(s::DeserializerState, ::Type{<:AffineScheme})
  R = load_typed_object(s, :coordinate_ring)
  return AffineScheme(R)
end

@register_serialization_type AffineVariety uses_id

function save_object(s::SerializerState, X::AffineVariety)
  save_data_dict(s) do
    save_typed_object(s, underlying_scheme(X), :underlying_scheme)
  end
end

function load_object(s::DeserializerState, ::Type{<:AffineVariety})
  X = load_typed_object(s, :underlying_scheme)
  return AffineVariety(X; check=false)
end

@register_serialization_type AffineAlgebraicSet uses_id

function save_object(s::SerializerState, X::AffineAlgebraicSet)
  save_data_dict(s) do
    save_typed_object(s, underlying_scheme(X), :underlying_scheme)
    # TODO: Also save the reduced version of X if any.
  end
end

function load_object(s::DeserializerState, ::Type{<:AffineAlgebraicSet})
  X = load_typed_object(s, :underlying_scheme)
  return AffineAlgebraicSet(X)
end

@register_serialization_type PrincipalOpenSubset uses_id

function save_object(s::SerializerState, X::PrincipalOpenSubset)
  save_data_dict(s) do
    save_typed_object(s, underlying_scheme(X), :underlying_scheme)
    save_typed_object(s, ambient_scheme(X), :ambient_scheme)
    save_typed_object(s, complement_equations(X), :comp_eqns)
  end
end

function load_object(s::DeserializerState, ::Type{<:PrincipalOpenSubset})
  U = load_typed_object(s, :underlying_scheme)
  eqns = load_typed_object(s, :comp_eqns)
  X = load_typed_object(s, :ambient_scheme)
  return PrincipalOpenSubset(X, U, eqns; check=false)
end

@register_serialization_type AffineSchemeMor uses_id

function save_object(s::SerializerState, f::AffineSchemeMor)
  save_data_dict(s) do
    save_typed_object(s, domain(f), :domain)
    save_typed_object(s, codomain(f), :codomain)
    save_typed_object(s, pullback(f), :pullback)
  end
end

function load_object(s::DeserializerState, ::Type{<:AffineSchemeMor})
  X = load_typed_object(s, :domain)
  Y = load_typed_object(s, :codomain)
  pbf = load_typed_object(s, :pullback)
  return AffineSchemeMor(X, Y, pbf; check=false)
end

@register_serialization_type SimpleGluing uses_id

function save_object(s::SerializerState, glue::SimpleGluing)
  save_data_dict(s) do
    A, B = patches(glue)
    f, g = gluing_morphisms(glue)
    save_typed_object(s, A, :first_patch)
    save_typed_object(s, B, :second_patch)
    save_typed_object(s, f, :first_map)
    save_typed_object(s, g, :second_map)
  end
end

function load_object(s::DeserializerState, ::Type{<:SimpleGluing})
  A = load_typed_object(s, :first_patch)
  B = load_typed_object(s, :second_patch)
  f = load_typed_object(s, :first_map)
  g = load_typed_object(s, :second_map)
  return SimpleGluing(A, B, f, g; check=false)
end

@register_serialization_type LazyGluing uses_id

function save_object(s::SerializerState, glue::LazyGluing)
  save_data_dict(s) do
    A, B = patches(glue)
    gd = glue.GD
    save_typed_object(s, A, :first_patch)
    save_typed_object(s, B, :second_patch)
    save_typed_object(s, gd, :gluing_data)
  end
end

function load_object(s::DeserializerState, ::Type{<:LazyGluing})
  A = load_typed_object(s, :first_patch)
  B = load_typed_object(s, :second_patch)
  gd = load_typed_object(s, :gluing_data)
  return LazyGluing(A, B, gd)
end

@register_serialization_type Covering uses_id

function save_object(s::SerializerState, cov::Covering)
  save_data_dict(s) do
    save_typed_object(s, patches(cov), :patches)
    # Sep. 18, 2024: Antony wants to refactor serialization of dictionaries.
    # When that is done, we should clean up this bit here.
    list = Tuple([(a, b) for (a, b) in gluings(cov)])
    save_typed_object(s, list, :gluings)
  end
end

function load_object(s::DeserializerState, ::Type{<:Covering})
  U = load_typed_object(s, :patches)
  list = load_typed_object(s, :gluings)
  dict = IdDict{Tuple{AbsAffineScheme, AbsAffineScheme}, AbsGluing}(a => b for (a, b) in list)
  return Covering(U, dict; check=false)
end

@register_serialization_type CoveredScheme uses_id

function save_object(s::SerializerState, X::CoveredScheme)
  save_data_dict(s) do
    save_typed_object(s, default_covering(X), :default_covering)
    # TODO: Look into whether or not to store refinements
    # It might be possible that we do not even need this, but we probably 
    # have to clarify the assumptions made throughout the code here.
  end
end

function load_object(s::DeserializerState, ::Type{<:CoveredScheme})
  cov = load_typed_object(s, :default_covering)
  return CoveredScheme(cov)
end

@register_serialization_type CoveringMorphism uses_id

function save_object(s::SerializerState, phi::CoveringMorphism)
  save_data_dict(s) do
    save_typed_object(s, domain(phi), :domain)
    save_typed_object(s, codomain(phi), :codomain)
    # TODO: Also clean up this part once the dicts are done.
    list = Tuple([(U, f) for (U, f) in morphisms(phi)])
    save_typed_object(s, list, :morphisms)
  end
end

function load_object(s::DeserializerState, ::Type{<:CoveringMorphism})
  X = load_typed_object(s, :domain)
  Y = load_typed_object(s, :codomain)
  list = load_typed_object(s, :morphisms)
  dict = IdDict{AbsAffineScheme, AbsAffineSchemeMor}(a => b for (a, b) in list)
  return CoveringMorphism(X, Y, dict; check=false)
end

@register_serialization_type CoveredSchemeMorphism uses_id

function save_object(s::SerializerState, phi::CoveredSchemeMorphism)
  save_data_dict(s) do
    save_typed_object(s, domain(phi), :domain)
    save_typed_object(s, codomain(phi), :codomain)
    save_typed_object(s, covering_morphism(phi), :covering_morphism)
  end
end

function load_object(s::DeserializerState, ::Type{<:CoveredSchemeMorphism})
  X = load_typed_object(s, :domain)
  Y = load_typed_object(s, :codomain)
  cov_mor = load_typed_object(s, :covering_morphism)
  return CoveredSchemeMorphism(X, Y, cov_mor)
end

