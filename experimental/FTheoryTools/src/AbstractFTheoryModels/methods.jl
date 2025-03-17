##########################################
### (1) Blow_up of F-theory models
##########################################


@doc raw"""
    blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")

Resolve an F-theory model by blowing up a locus in the ambient space.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> blow_up(t, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
Here is an example for a Weierstrass model.

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> blow_up(w, ["x", "y", "x1"]; coordinate_name = "e1")
Partially resolved Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)
```
"""
function blow_up(m::AbstractFTheoryModel, ideal_gens::Vector{String}; coordinate_name::String = "e")
  R = cox_ring(ambient_space(m))
  I = ideal([eval_poly(k, R) for k in ideal_gens])
  return blow_up(m, I; coordinate_name = coordinate_name)
end


@doc raw"""
    blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e")

Resolve an F-theory model by blowing up a locus in the ambient space.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4, x, y, z = gens(cox_ring(ambient_space(t)))
7-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4
 x
 y
 z

julia> blow_up(t, ideal([x, y, x1]); coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function blow_up(m::AbstractFTheoryModel, I::MPolyIdeal; coordinate_name::String = "e")
  return blow_up(m, ideal_sheaf(ambient_space(m), I); coordinate_name = coordinate_name)
end

function _ideal_sheaf_to_minimal_supercone_coordinates(X::AbsCoveredScheme, I::AbsIdealSheaf; coordinate_name::String = "e")
  # Return this when cannot convert ideal to minimal supercone coordinates
  not_possible = nothing

  # X needs to be a smooth toric variety
  X isa NormalToricVarietyType || return not_possible

  I isa ToricIdealSheafFromCoxRingIdeal || return not_possible
  defining_ideal = ideal_in_cox_ring(I)
  all(in(gens(base_ring(defining_ideal))), gens(defining_ideal)) || return not_possible
  R = cox_ring(X)
  coords = zeros(QQ, n_rays(X))
  for i in 1:n_rays(X)
    R[i] in gens(defining_ideal) && (coords[i] = 1)
  end
  coords == zeros(QQ, n_rays(X)) && return not_possible
  is_minimal_supercone_coordinate_vector(polyhedral_fan(X), coords) || return not_possible
  return coords
end

@doc raw"""
    blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e")

Resolve an F-theory model by blowing up a locus in the ambient space.
For this method, the blowup center is encoded by an ideal sheaf.

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> x1, x2, x3, x4, x, y, z = gens(cox_ring(ambient_space(t)))
7-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
 x1
 x2
 x3
 x4
 x
 y
 z

julia> blowup_center = ideal_sheaf(ambient_space(t), ideal([x, y, x1]))
Sheaf of ideals
  on normal, simplicial toric variety
with restrictions
   1: Ideal (x_5_1, x_4_1, x_1_1)
   2: Ideal (1)
   3: Ideal (x_5_3, x_4_3, x_1_3)
   4: Ideal (x_5_4, x_4_4, x_1_4)
   5: Ideal (1)
   6: Ideal (1)
   7: Ideal (1)
   8: Ideal (1)
   9: Ideal (1)
  10: Ideal (1)
  11: Ideal (1)
  12: Ideal (1)

julia> blow_up(t, blowup_center; coordinate_name = "e1")
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function blow_up(m::AbstractFTheoryModel, I::AbsIdealSheaf; coordinate_name::String = "e")
  
  # Cannot (yet) blowup if this is not a Tate or Weierstrass model
  entry_test = (m isa GlobalTateModel) || (m isa WeierstrassModel)
  @req entry_test "Blowups are currently only supported for Tate and Weierstrass models"
  @req (base_space(m) isa FamilyOfSpaces) == false "Base space must be a concrete space for blowups to work"

  # Compute the new ambient_space
  coords = _ideal_sheaf_to_minimal_supercone_coordinates(ambient_space(m), I)
  if !isnothing(coords)
    # Apply toric method
    bd = blow_up_along_minimal_supercone_coordinates(
      ambient_space(m), coords; coordinate_name=coordinate_name
    )
  else
    # Reroute to scheme theory
    bd = blow_up(I)
  end
  new_ambient_space = domain(bd)

  # Compute the new base
  # FIXME: THIS WILL IN GENERAL BE WRONG! IN PRINCIPLE, THE ABOVE ALLOWS TO BLOW UP THE BASE AND THE BASE ONLY.
  # FIXME: We should save the projection \pi from the ambient space to the base space.
  # FIXME: This is also ties in with the model sections to be saved, see below. Should the base change, so do these sections...
  new_base = base_space(m)

  # Construct the new model
  if m isa GlobalTateModel
    if isdefined(m, :tate_polynomial) && new_ambient_space isa NormalToricVariety
      f = tate_polynomial(m)
      new_tate_polynomial = strict_transform(bd, f)
      model = GlobalTateModel(explicit_model_sections(m), model_section_parametrization(m), new_tate_polynomial, base_space(m), new_ambient_space)
    else
      if bd isa ToricBlowupMorphism
        new_tate_ideal_sheaf = ideal_sheaf(domain(bd), strict_transform(bd, ideal_in_cox_ring(tate_ideal_sheaf(m))))
      else
        new_tate_ideal_sheaf = strict_transform(bd, tate_ideal_sheaf(m))
      end
      model = GlobalTateModel(explicit_model_sections(m), model_section_parametrization(m), new_tate_ideal_sheaf, base_space(m), new_ambient_space)
    end
  else
    if isdefined(m, :weierstrass_polynomial) && new_ambient_space isa NormalToricVariety
      f = weierstrass_polynomial(m)
      new_weierstrass_polynomial = strict_transform(bd, f)
      model = WeierstrassModel(explicit_model_sections(m), model_section_parametrization(m), new_weierstrass_polynomial, base_space(m), new_ambient_space)
    else
      if bd isa ToricBlowupMorphism
        new_weierstrass_ideal_sheaf = ideal_sheaf(domain(bd), strict_transform(bd, ideal_in_cox_ring(weierstrass_ideal_sheaf(m))))
      else
        new_weierstrass_ideal_sheaf = strict_transform(bd, weierstrass_ideal_sheaf(m))
      end
      model = WeierstrassModel(explicit_model_sections(m), model_section_parametrization(m), new_weierstrass_ideal_sheaf, base_space(m), new_ambient_space)
    end
  end

  # Copy/overwrite/set attributes
  model_attributes = m.__attrs
  for (key, value) in model_attributes
    set_attribute!(model, key, value)
  end
  set_attribute!(model, :partially_resolved, true)
  set_attribute!(model, :blow_down_morphism, bd)

  # Return the model
  return model
end



##########################################
### (2) Tuning
##########################################

# FIXME: The below tune function is not consistent with our other "tune" functions (it does not correspond mathematically to a tuning)
# and produces an undesired form for explicit_model_sections. To be improved at a later date.

# @doc raw"""
#     tune(m::AbstractFTheoryModel, p::MPolyRingElem; completeness_check::Bool = true)

# Tune an F-theory model by replacing the hypersurface equation by a custom (polynomial)
# equation. The latter can be any type of polynomial: a Tate polynomial, a Weierstrass
# polynomial or a general polynomial. We do not conduct checks to tell which type the
# provided polynomial is. Consequently, this tuning will always return a hypersurface model.

# Note that there is less functionality for hypersurface models than for Weierstrass or Tate
# models. For instance, `singular_loci` can (currently) not be computed for hypersurface models.

# # Examples
# ```jldoctest
# julia> B3 = projective_space(NormalToricVariety, 3)
# Normal toric variety

# julia> w = torusinvariant_prime_divisors(B3)[1]
# Torus-invariant, prime divisor on a normal toric variety

# julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
# Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

# Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

# julia> x1, x2, x3, x4, x, y, z = gens(parent(tate_polynomial(t)))
# 7-element Vector{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}:
#  x1
#  x2
#  x3
#  x4
#  x
#  y
#  z

# julia> new_tate_polynomial = x^3 - y^2 - x * y * z * x4^4
# -x4^4*x*y*z + x^3 - y^2

# julia> tuned_t = tune(t, new_tate_polynomial)
# Hypersurface model over a concrete base

# julia> hypersurface_equation(tuned_t) == new_tate_polynomial
# true

# julia> base_space(tuned_t) == base_space(t)
# true
# ```
# """
# function tune(m::AbstractFTheoryModel, p::MPolyRingElem; completeness_check::Bool = true)
#   entry_test = (m isa GlobalTateModel) || (m isa WeierstrassModel) || (m isa HypersurfaceModel)
#   @req entry_test "Tuning currently supported only for Weierstrass, Tate and hypersurface models"
#   @req (base_space(m) isa NormalToricVariety) "Currently, tuning is only supported for models over concrete toric bases"
#   equation = hypersurface_equation(m)
#   @req parent(p) == parent(equation) "Parent mismatch between given and existing hypersurface polynomial"
#   @req degree(p) == degree(equation) "Degree mismatch between given and existing hypersurface polynomial"
#   p == equation && return m
#   explicit_model_sections = Dict{String, MPolyRingElem}()
#   gens_S = gens(parent(p))
#   for k in 1:length(gens_S)
#     explicit_model_sections[string(gens_S[k])] = gens_S[k]
#   end
#   tuned_model = HypersurfaceModel(explicit_model_sections, p, p, base_space(m), ambient_space(m), fiber_ambient_space(m))
#   set_attribute!(tuned_model, :partially_resolved, false)
#   return tuned_model
# end



@doc raw"""
    put_over_concrete_base(m::AbstractFTheoryModel, concrete_data::Dict{String, <:Any}; completeness_check::Bool = true)

Put an F-theory model defined over a family of spaces over a concrete base.

Currently, this functionality is limited to Tate and Weierstrass models.

# Examples
```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", completeness_check = false)
Assuming that the first row of the given grading is the grading under Kbar

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
        @req parent(concrete_data[gen_name]) == cox_ring(concrete_data["base"]) "Specified sections must reside in Cox ring of given base"
        new_model_secs[gen_name] = concrete_data[gen_name]
      end
    end

    # Create ring map to evaluate Weierstrass/Tate sections
    parametrization = model_section_parametrization(m)
    domain = parent(collect(values(parametrization))[1])
    codomain = cox_ring(concrete_data["base"])
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



##########################################
### (3) Meta data setters
##########################################

function set_arxiv_id(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_id => desired_value)
end

function set_arxiv_doi(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_doi => desired_value)
end

function set_arxiv_link(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_link => desired_value)
end

function set_arxiv_model_equation_number(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_model_equation_number => desired_value)
end

function set_arxiv_model_page(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_model_page => desired_value)
end

function set_arxiv_model_section(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_model_section => desired_value)
end

function set_arxiv_version(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :arxiv_version => desired_value)
end

function set_associated_literature_models(m::AbstractFTheoryModel, desired_value::Vector{String})
  set_attribute!(m, :associated_literature_models => desired_value)
end

function set_journal_doi(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_doi => desired_value)
end

function set_journal_link(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_link => desired_value)
end

function set_journal_model_equation_number(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_model_equation_number => desired_value)
end

function set_journal_model_page(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_model_page => desired_value)
end

function set_journal_model_section(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_model_section => desired_value)
end

function set_journal_name(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_name => desired_value)
end

function set_journal_pages(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_pages => desired_value)
end

function set_journal_report_numbers(m::AbstractFTheoryModel, desired_value::Vector{String})
  set_attribute!(m, :journal_report_numbers => desired_value)
end

function set_journal_volume(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_volume => desired_value)
end

function set_journal_year(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :journal_year => desired_value)
end

function set_literature_identifier(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :literature_identifier => desired_value)
end

function set_model_description(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :model_description => desired_value)
end

function set_model_parameters(m::AbstractFTheoryModel, desired_value::Vector{String})
  set_attribute!(m, :model_parameters => desired_value)
end

function set_paper_authors(m::AbstractFTheoryModel, desired_value::Vector{String})
  set_attribute!(m, :paper_authors => desired_value)
end

function set_paper_buzzwords(m::AbstractFTheoryModel, desired_value::Vector{String})
  set_attribute!(m, :paper_buzzwords => desired_value)
end

function set_paper_description(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :paper_description => desired_value)
end

function set_paper_title(m::AbstractFTheoryModel, desired_value::String)
  set_attribute!(m, :paper_title => desired_value)
end

function set_birational_literature_models(m::AbstractFTheoryModel, desired_value::Vector{String})
  set_attribute!(m, :birational_literature_models => desired_value)
end



##########################################
### (4) Meta data adders
##########################################

function add_associated_literature_model(m::AbstractFTheoryModel, addition::String)
  values = has_associated_literature_models(m) ? associated_literature_models(m) : []
  !(addition in values) && set_attribute!(m, :associated_literature_models => vcat(values, [addition]))
end

function add_journal_report_number(m::AbstractFTheoryModel, addition::String)
  values = has_journal_report_numbers(m) ? journal_report_numbers(m) : []
  !(addition in values) && set_attribute!(m, :journal_report_numbers => vcat(values, [addition]))
end

function add_model_parameter(m::AbstractFTheoryModel, addition::String)
  values = has_model_parameters(m) ? model_parameters(m) : []
  !(addition in values) && set_attribute!(m, :model_parameters => vcat(values, [addition]))
end

function add_paper_author(m::AbstractFTheoryModel, addition::String)
  values = has_paper_authors(m) ? paper_authors(m) : []
  !(addition in values) && set_attribute!(m, :paper_authors => vcat(values, [addition]))
end

function add_paper_buzzword(m::AbstractFTheoryModel, addition::String)
  values = has_paper_buzzwords(m) ? paper_buzzwords(m) : []
  !(addition in values) && set_attribute!(m, :paper_buzzwords => vcat(values, [addition]))
end

function add_birational_literature_model(m::AbstractFTheoryModel, addition::String)
  values = has_birational_literature_models(m) ? birational_literature_models(m) : []
  !(addition in values) && set_attribute!(m, :birational_literature_models => vcat(values, [addition]))
end



##########################################
### (5) Specialized model data setters
##########################################

function set_generating_sections(m::AbstractFTheoryModel, vs::Vector{Vector{String}})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  set_attribute!(m, :generating_sections => [[f(eval_poly(l, R)) for l in k] for k in vs])
end

function set_torsion_sections(m::AbstractFTheoryModel, vs::Vector{Vector{String}})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  set_attribute!(m, :torsion_sections => [[f(eval_poly(l, R)) for l in k] for k in vs])
end

function set_resolutions(m::AbstractFTheoryModel, desired_value::Vector{Tuple{Vector{Vector{String}}, Vector{String}}})
  set_attribute!(m, :resolutions => desired_value)
end

function set_resolution_generating_sections(m::AbstractFTheoryModel, vs::Vector{Vector{Vector{Vector{String}}}})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  result = [[[[f(eval_poly(a, R)) for a in b] for b in c] for c in d] for d in vs]
  set_attribute!(m, :resolution_generating_sections => result)
end

function set_resolution_zero_sections(m::AbstractFTheoryModel, vs::Vector{Vector{Vector{String}}})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  result = [[[f(eval_poly(a, R)) for a in b] for b in c] for c in vs]
  set_attribute!(m, :resolution_zero_sections => result)
end

function set_weighted_resolutions(m::AbstractFTheoryModel, desired_value::Vector{Tuple{Vector{Tuple{Vector{String}, Vector{Int}}}, Vector{String}}})
  set_attribute!(m, :weighted_resolutions => desired_value)
end

function set_weighted_resolution_generating_sections(m::AbstractFTheoryModel, vs::Vector{Vector{Vector{Vector{String}}}})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  result = [[[[f(eval_poly(a, R)) for a in b] for b in c] for c in d] for d in vs]
  set_attribute!(m, :weighted_resolution_generating_sections => result)
end

function set_weighted_resolution_zero_sections(m::AbstractFTheoryModel, vs::Vector{Vector{Vector{String}}})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  result = [[[f(eval_poly(a, R)) for a in b] for b in c] for c in vs]
  set_attribute!(m, :weighted_resolution_zero_sections => result)
end

function set_zero_section(m::AbstractFTheoryModel, desired_value::Vector{String})
  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(m))), cached = false)
  f = hom(R, cox_ring(base_space(m)), collect(values(explicit_model_sections(m))))
  set_attribute!(m, :zero_section => [f(eval_poly(l, R)) for l in desired_value])
end

function set_zero_section_class(m::AbstractFTheoryModel, desired_value::String)
  divs = torusinvariant_prime_divisors(ambient_space(m))
  cohomology_ring(ambient_space(m); check=false)
  cox_gens = string.(gens(cox_ring(ambient_space(m))))
  @req desired_value in cox_gens "Specified zero section is invalid"
  index = findfirst(==(desired_value), cox_gens)
  set_attribute!(m, :zero_section_index => index::Int)
  set_attribute!(m, :zero_section_class => cohomology_class(divs[index]))
end

function set_gauge_algebra(m::AbstractFTheoryModel, algebras::Vector{String})
  C = algebraic_closure(QQ)
  function _construct(g::String)
    if g == "0"
      return abelian_lie_algebra(C, 0)
    elseif g == "u(1)"
      return lie_algebra(C,1,[C(1im)*identity_matrix(C,1)],["i"])
    elseif g[1:2] == "su"
      return special_linear_lie_algebra(C, parse(Int, g[4:end-1]))
    elseif g[1:2] == "so"
      return special_orthogonal_lie_algebra(C, parse(Int, g[4:end-1]))
    elseif g[1:2] == "sp"
      return symplectic_lie_algebra(C, parse(Int, g[4:end-1]))
    elseif g[1:1] == "e"
      return lie_algebra(C, :E, parse(Int, g[3:end-1]))
    elseif g[1:1] == "f"
      return lie_algebra(C, :F, parse(Int, g[3:end-1]))
    elseif g[1:1] == "g"
      return lie_algebra(C, :G, parse(Int, g[3:end-1]))
    end
    error("Unknown algebra description")
  end
  set_attribute!(m, :gauge_algebra => direct_sum(C, LieAlgebra{elem_type(C)}[_construct(g) for g in algebras]))
end

function set_global_gauge_quotients(m::AbstractFTheoryModel, quotients::Vector{Vector{String}})
 set_attribute!(m, :global_gauge_quotients => quotients)
end


##########################################
### (6) Specialized model data adders
##########################################

function add_generating_section(m::AbstractFTheoryModel, addition::Vector{String})
  values = has_generating_sections(m) ? birational_literature_models(m) : []
  !(addition in values) && set_attribute!(m, :generating_sections => vcat(values, [addition]))
end

@doc raw"""
    add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})

Add a known resolution for a model.

```jldoctest
julia> m = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> add_resolution(m, [["x", "y"], ["y", "s", "w"], ["s", "e4"], ["s", "e3"], ["s", "e1"]], ["s", "w", "e3", "e1", "e2"])

julia> length(resolutions(m))
2
```
"""
function add_resolution(m::AbstractFTheoryModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
  @req length(exceptionals) == length(centers) "Number of exceptionals must match number of centers"
  resolution = [centers, exceptionals]
  known_resolutions = has_resolutions(m) ? resolutions(m) : []
  !(resolution in known_resolutions) && set_attribute!(m, :resolutions => vcat(known_resolutions, [resolution]))
end

function add_resolution_generating_section(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
  values = has_resolution_generating_sections(m) ? resolution_generating_sections(m) : []
  !(addition in values) && set_attribute!(m, :resolution_generating_sections => vcat(values, [addition]))
end

function add_resolution_zero_section(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
  values = has_resolution_zero_sections(m) ? resolution_zero_sections(m) : []
  !(addition in values) && set_attribute!(m, :resolution_zero_sections => vcat(values, [addition]))
end

function add_weighted_resolution(m::AbstractFTheoryModel, addition::Vector{Vector})
  values = has_weighted_resolutions(m) ? weighted_resolutions(m) : []
  !(addition in values) && set_attribute!(m, :weighted_resolutions => vcat(values, [addition]))
end

function add_weighted_resolution_generating_section(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
  values = has_weighted_resolution_generating_sections(m) ? weighted_resolution_generating_sections(m) : []
  !(addition in values) && set_attribute!(m, :weighted_resolution_generating_sections => vcat(values, [addition]))
end

function add_weighted_resolution_zero_section(m::AbstractFTheoryModel, addition::Vector{Vector{Vector{String}}})
  values = has_weighted_resolution_zero_sections(m) ? weighted_resolution_zero_sections(m) : []
  !(addition in values) && set_attribute!(m, :weighted_resolution_zero_sections => vcat(values, [addition]))
end



##########################################
### (7) Specialized model methods
##########################################

@doc raw"""
    resolve(m::AbstractFTheoryModel, index::Int)

Resolve a model with the index-th resolution that is known.

```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> t2 = resolve(t, 1)
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> cox_ring(ambient_space(t2))
Multivariate polynomial ring in 12 variables over QQ graded by
  x1 -> [1 0 0 0 0 0 0]
  x2 -> [0 1 0 0 0 0 0]
  x3 -> [0 1 0 0 0 0 0]
  x4 -> [0 1 0 0 0 0 0]
  x -> [0 0 1 0 0 0 0]
  y -> [0 0 0 1 0 0 0]
  z -> [0 0 0 0 1 0 0]
  e1 -> [0 0 0 0 0 1 0]
  e4 -> [0 0 0 0 0 0 1]
  e2 -> [-1 -3 -1 1 -1 -1 0]
  e3 -> [0 4 1 -1 1 0 -1]
  s -> [2 6 -1 0 2 1 1]

julia> w2 = 2 * torusinvariant_prime_divisors(B3)[1]
Torus-invariant, non-prime divisor on a normal toric variety

julia> t3 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w2), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> t4 = resolve(t3, 1)
Partially resolved global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```
"""
function resolve(m::AbstractFTheoryModel, resolution_index::Int)

  # For model 1511.03209 and resolution_index = 1, a particular resolution is available from an artifact
  if has_attribute(m, :arxiv_id)
    if resolution_index == 1 && arxiv_id(m) == "1511.03209"

      # Once zenodo is updated, the following two lines can be executed to cut the hour taking resolution procedure down to minutes
      #=
      model_data_path = artifact"FTM-1511-03209/1511-03209-resolved.mrdi"
      return load(model_data_path)
      =#

    end
  end

  # To be extended to hypersurface models...
  entry_test = (m isa GlobalTateModel) || (m isa WeierstrassModel)
  @req entry_test "Resolve currently supported only for Weierstrass and Tate models"
  @req (base_space(m) isa NormalToricVariety) "Currently, resolve is only supported for models over concrete toric bases"
  @req (ambient_space(m) isa NormalToricVariety) "Currently, resolve is only supported for singular models defined in a toric space"
  @req has_attribute(m, :resolutions) "No resolutions known for this model"
  @req resolution_index > 0 "The resolution must be specified by a non-negative integer"
  @req resolution_index <= length(resolutions(m)) "The resolution must be specified by an integer that is not larger than the number of known resolutions"
  
  # Gather information for resolution
  centers, exceptionals = resolutions(m)[resolution_index]
  nr_blowups = length(centers)
  
  # Resolve the model
  resolved_model = m
  blow_up_chain = []
  for k in 1:nr_blowups

    # Replace parameters in the blow_up_center with explicit_model_sections
    blow_up_center = centers[k]
    for l in 1:length(blow_up_center)
      model_sections = explicit_model_sections(resolved_model)
      if haskey(model_sections, blow_up_center[l])
        new_locus = string(explicit_model_sections(resolved_model)[blow_up_center[l]])
        blow_up_center[l] = new_locus
      end
    end

    # Conduct the blowup
    if ambient_space(resolved_model) isa NormalToricVariety
      # Toric case is easy...
      resolved_model = blow_up(resolved_model, blow_up_center; coordinate_name = exceptionals[k])
    else
      
      # Compute proper transform of center generated by anything but exceptional divisors
      filtered_center = [c for c in blow_up_center if !(c in exceptionals)]
      initial_ambient_space = ambient_space(m)
      initial_cox_ring = cox_ring(initial_ambient_space)
      initial_filtered_ideal_sheaf = ideal_sheaf(initial_ambient_space, ideal([eval_poly(l, initial_cox_ring) for l in filtered_center]))
      bd_morphism = get_attribute(blow_up_chain[1], :blow_down_morphism)
      filtered_ideal_sheaf = strict_transform(bd_morphism, initial_filtered_ideal_sheaf)
      for l in 2:k-1
        bd_morphism = get_attribute(blow_up_chain[l], :blow_down_morphism)
        filtered_ideal_sheaf = strict_transform(bd_morphism, filtered_ideal_sheaf)
      end

      # Compute strict transform of ideal sheaves appearing in blowup center
      exceptional_center = [c for c in blow_up_center if (c in exceptionals)]
      positions = [findfirst(==(l), exceptionals) for l in exceptional_center]
      exceptional_divisors = [exceptional_divisor(get_attribute(blow_up_chain[l], :blow_down_morphism)) for l in positions]
      exceptional_ideal_sheafs = [ideal_sheaf(d) for d in exceptional_divisors]
      for l in 1:length(positions)
        if positions[l] < k-1
          for m in positions[l]+1: k-1
            internal_bd_morphism = get_attribute(blow_up_chain[m], :blow_down_morphism)
            exceptional_ideal_sheafs[l] = strict_transform(internal_bd_morphism, exceptional_ideal_sheafs[l])
          end
        end
      end

      # Compute the prepared center
      prepared_center = filtered_ideal_sheaf
      if length(exceptional_ideal_sheafs) > 0
        prepared_center = prepared_center + sum(exceptional_ideal_sheafs)
      end

      # Execute the blow-up
      resolved_model = blow_up(resolved_model, prepared_center; coordinate_name = exceptionals[k])
    end

    # Remember the result
    push!(blow_up_chain, resolved_model)

  end


  # For model 1511.03209 and resolution_index = 1, we extend beyond what is currently saved as resolution in our json file.
  # Once zenodo artifact is updated, comment out the following lines
  if has_attribute(m, :arxiv_id)
    if resolution_index == 1 && arxiv_id(m) == "1511.03209"

      # Bring in previous result. For time reasons, we can also load this in. This line with the hard-coded path is to be removed eventually...
      resolved_model = load("/datadisk/ToDo/1511.03209/Update/new_1511.03209-resolved.mrdi")
      
      # Now execute 3 more blowups

      # Additional blowup 1:
      # Additional blowup 1:

      # Ambient space has the following rays
      #x: [0,0,0,-3,1]
      #y: [0,0,0,2,-1]
      #z: [0,0,0,0,1]
      # We add the ray m1: (0,0,0,1,0). This looks like y + z = 2 * m1.
      # So naively, I think of this as blowing up y^2 = z^2 = 0 and introducing the variable m1. For the strict transform, we thus do
      # y^2 -> y^2 * m1
      # z^2 -> z^2 * m1
      # y * z -> y * z * m1
      as = ambient_space(resolved_model);
      bl = domain(blow_up(as, [0,0,0,1,0], coordinate_name = "m1"));
      f = hypersurface_equation(resolved_model);
      my_mons = collect(monomials(f));
      pos_1 = findfirst(k -> k == "y", [string(a) for a in gens(cox_ring(as))])
      pos_2 = findfirst(k -> k == "z", [string(a) for a in gens(cox_ring(as))])
      exp_list = [collect(exponents(m))[1] for m in my_mons];
      my_exps = [[k[pos_1], k[pos_2]] for k in exp_list];
      @req all(k -> isinteger(1//2 * sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
      m_power = [Int(1//2 * sum(a)) for a in my_exps]
      overall_factor = minimum(m_power)
      new_coeffs = collect(coefficients(f))
      new_exps = [vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)]
      my_builder = MPolyBuildCtx(cox_ring(bl))
      for a in 1:length(new_exps)
        push_term!(my_builder, new_coeffs[a], new_exps[a])
      end
      new_tate_polynomial = finish(my_builder);
      model_bl = GlobalTateModel(explicit_model_sections(resolved_model), model_section_parametrization(resolved_model), new_tate_polynomial, base_space(resolved_model), bl);
      set_attribute!(model_bl, :partially_resolved, true)
      
      
      # Additional blowup 2:
      # Additional blowup 2:

      # Ambient space has the following rays:
      # x: [0,0,0,-3,1]
      # y: [0,0,0,2,-1]
      # z: [0,0,0,0,1]
      # m1: [0,0,0,1,0]
      # We add the ray m2: (0,0,0,-2,0). This looks like blowing up at x = m1 = 0 and introducing the variable m2. For the strict transform, we thus do
      # x -> x * m2
      # m1 -> m1 * m2
      as = ambient_space(model_bl);
      bl = domain(blow_up(as, [0,0,0,-2,1], coordinate_name = "m2"));
      f = hypersurface_equation(model_bl);
      my_mons = collect(monomials(f));
      pos_1 = findfirst(k -> k == "x", [string(a) for a in gens(cox_ring(as))])
      pos_2 = findfirst(k -> k == "m1", [string(a) for a in gens(cox_ring(as))])
      exp_list = [collect(exponents(m))[1] for m in my_mons];
      my_exps = [[k[pos_1], k[pos_2]] for k in exp_list];
      @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
      m_power = [Int(sum(a)) for a in my_exps]
      overall_factor = minimum(m_power)
      new_coeffs = collect(coefficients(f))
      new_exps = [vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)]
      my_builder = MPolyBuildCtx(cox_ring(bl))
      for a in 1:length(new_exps)
        push_term!(my_builder, new_coeffs[a], new_exps[a])
      end
      new_tate_polynomial = finish(my_builder);
      model_bl2 = GlobalTateModel(explicit_model_sections(model_bl), model_section_parametrization(model_bl), new_tate_polynomial, base_space(model_bl), bl);
      set_attribute!(model_bl2, :partially_resolved, true)


      # Additional blowup 3:
      # Additional blowup 3:
      
      # Ambient space has the following rays:
      # x: [0,0,0,-3,1]
      # y: [0,0,0,2,-1]
      # z: [0,0,0,0,1]
      # m1: [0,0,0,1,0]
      # m2: [0,0,0,-2,1]
      # We add the ray m3: (0,0,0,-1,1). This looks like blowing up at m1 = m2 = 0 and introducing the variable m3. For the strict transform, we thus do
      # m1 -> m1 * m3
      # m2 -> m2 * m3
      # (0,0,0,-1,1) = (0,0,0,-3,1) + 2 * (0,0,0,1,0), so also x * m1^2 * m3 is homogeneous.

      as = ambient_space(model_bl2);

      # Q: Which maximal cones of as do contain the ray (0,0,0,-1,1)? For those are then telling us which variables need scaling.
      # -> This tells us the overall pwer of m3!
      Sigma = polyhedral_fan(as)
      Oscar.pm_object(Sigma).FACET_NORMALS
      Oscar.pm_object(Sigma).MAXIMAL_CONES_FACETS
      my_list = findall(mc -> [0,0,0,-1,1] in mc, maximal_cones(Sigma))
      # This vector appears in 198 maximal cones! Puh... how to take all of this into account?


      bl = domain(blow_up(as, [0,0,0,-1,1], coordinate_name = "m3"));
      f = hypersurface_equation(model_bl2);
      my_mons = collect(monomials(f));
      pos_1 = findfirst(k -> k == "m1", [string(a) for a in gens(cox_ring(as))])
      pos_2 = findfirst(k -> k == "m2", [string(a) for a in gens(cox_ring(as))])
      exp_list = [collect(exponents(m))[1] for m in my_mons];
      my_exps = [[k[pos_1], k[pos_2]] for k in exp_list];
      @req all(k -> isinteger(sum(k)), my_exps) "Inconsistency encountered in computation of strict transform. Please inform the authors."
      m_power = [Int(sum(a)) for a in my_exps]
      overall_factor = minimum(m_power)
      new_coeffs = collect(coefficients(f))
      new_exps = [vcat([exp_list[k], m_power[k] - overall_factor]...) for k in 1:length(exp_list)]
      #=
      my_builder = MPolyBuildCtx(cox_ring(bl))
      for a in 1:length(new_exps)
        push_term!(my_builder, new_coeffs[a], new_exps[a])
      end
      new_tate_polynomial = finish(my_builder);
      =#

      my_builder = MPolyBuildCtx(cox_ring(bl))
      for a in 1:Int(floor(1//2*length(new_exps)))
        push_term!(my_builder, new_coeffs[a], new_exps[a])
      end
      new_tate_polynomial_1 = finish(my_builder);

      my_builder = MPolyBuildCtx(cox_ring(bl))
      for a in Int(floor(1//2*length(new_exps)))+1:length(new_exps)
        push_term!(my_builder, new_coeffs[a], new_exps[a])
      end
      new_tate_polynomial_2 = finish(my_builder);
      @req is_homogeneous(new_tate_polynomial_1) "First poly not homogeneous"
      @req is_homogeneous(new_tate_polynomial_2) "Second poly not homogeneous"
      @req degree(new_tate_polynomial_2) == degree(new_tate_polynomial_1) "Poly have different degree"

      model_bl3 = GlobalTateModel(explicit_model_sections(model_bl2), model_section_parametrization(model_bl2), new_tate_polynomial, base_space(model_bl2), bl);
      set_attribute!(model_bl3, :partially_resolved, true)


      # We confirm that after these steps, we achieve what we desire.
      @req is_smooth(ambient_space(model_bl3)) "Ambient space not yet smooth. Damn it!"
      @req is_homogeneous(hypersurface_equation(model_bl3)) "Strict transform is not homogeneous. Damn it!"
      
      
      # Finally return the result
      return model_bl3

    end
  end


  return resolved_model
end
