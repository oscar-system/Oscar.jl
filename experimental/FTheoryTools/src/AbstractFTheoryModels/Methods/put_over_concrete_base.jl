@doc raw"""
    put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any}; completeness_check::Bool = true)

Put an F-theory model defined over a family of spaces over a concrete base.

Currently, this functionality is limited to Tate and Weierstrass models.

# Examples
```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", completeness_check = false)
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

julia> t2 = put_over_concrete_base(t, Dict("base" => B3, "w" => w, "a21" => a21, "a32" => a32, "a43" => a43), completeness_check = false)
Global Tate model over a concrete base
```
"""
function put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any}; completeness_check::Bool = true)
  # Conduct elementary entry checks
  @req base_space(m) isa FamilyOfSpaces "The model must be defined over a family of spaces"
  @req haskey(concrete_data, "base") "The base space must be specified"
  @req (concrete_data["base"] isa NormalToricVariety) "Currently, models over families of spaces can only be put over toric bases"
  @req ((m isa WeierstrassModel) || (m isa GlobalTateModel)) "Currently, only Tate or Weierstrass models can be put on a concrete base"
  
  # Work out the Weierstrass/Tate sections
  new_model_secs = Dict{String, MPolyRingElem}()
  if is_empty(model_section_parametrization(m))
    
    # No parametrization, so simply take generic sections
    
    # Pick generic values
    if m isa WeierstrassModel
      new_model_secs["f"] = generic_section(anticanonical_bundle(concrete_data["base"])^4)
      new_model_secs["g"] = generic_section(anticanonical_bundle(concrete_data["base"])^6)
    else
      new_model_secs["a1"] = generic_section(anticanonical_bundle(concrete_data["base"]))
      new_model_secs["a2"] = generic_section(anticanonical_bundle(concrete_data["base"])^2)
      new_model_secs["a3"] = generic_section(anticanonical_bundle(concrete_data["base"])^3)
      new_model_secs["a4"] = generic_section(anticanonical_bundle(concrete_data["base"])^4)
      new_model_secs["a6"] = generic_section(anticanonical_bundle(concrete_data["base"])^6)
    end
    
  else

    # Parametrization for Weierstrass/Tate sections found

    # Have all parametrizing sections been provided by the user?
    polys = collect(values(model_section_parametrization(m)))
    all_appearing_monomials = vcat([collect(monomials(p)) for p in polys]...)
    all_appearing_exponents = hcat([collect(exponents(m))[1] for m in all_appearing_monomials]...)
    for k in 1:nrows(all_appearing_exponents)
      if any(!is_zero, all_appearing_exponents[k,:])
        gen_name = string(symbols(parent(polys[1]))[k])
        @req haskey(concrete_data, gen_name) "Required base section $gen_name not specified"
        @req parent(concrete_data[gen_name]) == coordinate_ring(concrete_data["base"]) "Specified sections must reside in Cox ring of given base"
        new_model_secs[gen_name] = concrete_data[gen_name]
      end
    end

    # Create ring map to evaluate Weierstrass/Tate sections
    parametrization = model_section_parametrization(m)
    domain = parent(collect(values(parametrization))[1])
    codomain = coordinate_ring(concrete_data["base"])
    images = [haskey(new_model_secs, string(k)) ? new_model_secs[string(k)] : zero(codomain) for k in gens(domain)]
    mapper = hom(domain, codomain, images)

    # Compute defining sections
    if m isa WeierstrassModel

      if haskey(parametrization, "f")
        new_sec = mapper(parametrization["f"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"])^4)) "Degree mismatch"
        end
        new_model_secs["f"] = new_sec
      else
        new_model_secs["f"] = generic_section(anticanonical_bundle(concrete_data["base"])^4)
      end

      if haskey(parametrization, "g")
        new_sec = mapper(parametrization["g"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"])^6)) "Degree mismatch"
        end
        new_model_secs["g"] = new_sec
      else
        new_model_secs["g"] = generic_section(anticanonical_bundle(concrete_data["base"])^6)
      end

    else

      if haskey(parametrization, "a1")
        new_sec = mapper(parametrization["a1"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"]))) "Degree mismatch"
        end
        new_model_secs["a1"] = new_sec
      else
        new_model_secs["a1"] = generic_section(anticanonical_bundle(concrete_data["base"]))
      end

      if haskey(parametrization, "a2")
        new_sec = mapper(parametrization["a2"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"])^2)) "Degree mismatch"
        end
        new_model_secs["a2"] = new_sec
      else
        new_model_secs["a2"] = generic_section(anticanonical_bundle(concrete_data["base"])^2)
      end

      if haskey(parametrization, "a3")
        new_sec = mapper(parametrization["a3"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"])^3)) "Degree mismatch"
        end
        new_model_secs["a3"] = new_sec
      else
        new_model_secs["a3"] = generic_section(anticanonical_bundle(concrete_data["base"])^3)
      end

      if haskey(parametrization, "a4")
        new_sec = mapper(parametrization["a4"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"])^4)) "Degree mismatch"
        end
        new_model_secs["a4"] = new_sec
      else
        new_model_secs["a4"] = generic_section(anticanonical_bundle(concrete_data["base"])^4)
      end

      if haskey(parametrization, "a6")
        new_sec = mapper(parametrization["a6"])
        if !is_zero(new_sec)
          @req degree(new_sec) == degree(generic_section(anticanonical_bundle(concrete_data["base"])^6)) "Degree mismatch"
        end
        new_model_secs["a6"] = new_sec
      else
        new_model_secs["a6"] = generic_section(anticanonical_bundle(concrete_data["base"])^6)
      end

    end
  
  end

  # Compute the new model
  if m isa WeierstrassModel
    return weierstrass_model(concrete_data["base"], new_model_secs, model_section_parametrization(m); completeness_check)
  else
    return global_tate_model(concrete_data["base"], new_model_secs, model_section_parametrization(m); completeness_check)
  end
end
