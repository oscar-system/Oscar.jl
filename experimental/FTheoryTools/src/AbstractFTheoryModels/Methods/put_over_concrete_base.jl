@doc raw"""
    put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any})

Put an F-theory model defined over a family of spaces over a concrete base.

Currently, this functionality is limited to Tate and Weierstrass models.

Internally, this function uses generic sections of line bundles.

!!! note "Randomness"
    The random source used for randomized computations can be set with the `rng` keyword.

!!! note "Complete toric base"
    This function assumes that the toric base space is **complete**.
    Checking completeness may take a long time. To skip this check,
    pass the **optional keyword argument** `completeness_check=false`.

# Examples
```jldoctest
julia> using Random;

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w_bundle = toric_line_bundle(torusinvariant_prime_divisors(B3)[1])
Toric line bundle on a normal toric variety

julia> kbar = anticanonical_bundle(B3)
Toric line bundle on a normal toric variety

julia> w = generic_section(w_bundle);

julia> a21 = generic_section(kbar^2 * w_bundle^(-1));

julia> a32 = generic_section(kbar^3 * w_bundle^(-2));

julia> a43 = generic_section(kbar^4 * w_bundle^(-3));

julia> t2 = put_over_concrete_base(t, Dict("base" => B3, "w" => w, "a21" => a21, "a32" => a32, "a43" => a43), completeness_check = false, rng = Random.Xoshiro(1234))
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function put_over_concrete_base(
  m::AbstractFTheoryModel,
  concrete_data::ConcreteModelDataType;
  completeness_check::Bool=true,
  rng::AbstractRNG=Random.default_rng(),
)
  # Conduct elementary entry checks
  @req base_space(m) isa FamilyOfSpaces "The model must be defined over a family of spaces"
  @req haskey(concrete_data, "base") "The base space must be specified"
  @req (
    concrete_data["base"] isa NormalToricVariety
  ) "Currently, models over families of spaces can only be put over toric bases"
  @req (
    m isa Union{WeierstrassModel,GlobalTateModel}
  ) "Currently, only Tate or Weierstrass models can be put on a concrete base"

  base = concrete_data["base"]::NormalToricVariety
  parametrization = model_section_parametrization(m)
  defining_sections = if m isa WeierstrassModel
    (("f", 4), ("g", 6))
  else
    (("a1", 1), ("a2", 2), ("a3", 3), ("a4", 4), ("a6", 6))
  end
  anticanonical = anticanonical_bundle(base)
  new_model_sections = Dict{String,MPolyRingElem}()

  if isempty(parametrization)
    for (name, power) in defining_sections
      new_model_sections[name] = generic_section(anticanonical^power; rng)
    end
  else
    parametrization_ring = parent(first(values(parametrization)))
    parametrization_symbols = symbols(parametrization_ring)
    used_generators = falses(ngens(parametrization_ring))
    for polynomial in values(parametrization)
      for index in eachindex(used_generators)
        used_generators[index] |= degree(polynomial, index) > 0
      end
    end

    for index in eachindex(used_generators)
      used_generators[index] || continue
      name = string(parametrization_symbols[index])
      @req haskey(concrete_data, name) "Required base section " * name * " not specified"
      section = concrete_data[name]
      @req (
        parent(section) == coordinate_ring(base)
      ) "Specified sections must reside in Cox ring of given base"
      new_model_sections[name] = section
    end

    base_ring = coordinate_ring(base)
    zero_section = zero(base_ring)
    images = [
      get(new_model_sections, string(generator), zero_section) for
      generator in gens(parametrization_ring)
    ]
    evaluation_map = hom(parametrization_ring, base_ring, images)

    for (name, power) in defining_sections
      bundle = anticanonical^power
      if haskey(parametrization, name)
        section = evaluation_map(parametrization[name])
        if !is_zero(section)
          expected_degree = divisor_class(toric_divisor_class(bundle))
          @req degree(section) == expected_degree "Degree mismatch"
        end
        new_model_sections[name] = section
      else
        new_model_sections[name] = generic_section(bundle; rng)
      end
    end
  end

  result = if m isa WeierstrassModel
    weierstrass_model(
      base, new_model_sections, parametrization; completeness_check
    )
  else
    global_tate_model(
      base, new_model_sections, parametrization; completeness_check
    )
  end
  _copy_model_metadata!(result, m)
  return result
end
