#######################################################
# 1. User interface for literature models
#######################################################

@doc raw"""
    literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)

Many models have been created in the F-theory literature.
A significant number of them have even been given specific
names, for instance the "U(1)-restricted SU(5)-GUT model".
This method has access to a database, from which it can
look up such literature models.

Currently, you can provide
any combination of the following optional arguments
to the method `literature_model`:
* `doi`: A string representing the DOI of the publication that
introduced the model in question.
* `equation`: A string representing the number of the equation that introduced
the model in question.
For papers, that were posted on the arXiv, we can instead of the `doi` also
provide the following:
* `arxiv_id`: A string that represents the arXiv identifier of the paper that
introduced the model in question.
* `version`: A string representing the version of the arXiv upload.
The method `literature_model` attempts to find a model in our database
for which the provided data matches the information in our record. If no such
model can be found, or multiple models exist with information matching the
provided information, then an error is raised.

Some literature models require additional parameters to specified to single out
a model from a family of models. Such models can be provided using the optional
argument `model_parameters`, which should be a dictionary such as `Dict("k" => 5)`.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> v = ambient_space(t)
A family of spaces of dimension d = 5

julia> coordinate_ring(v)
Multivariate polynomial ring in 8 variables w, a1, a21, a32, ..., z
  over rational field
```
It is also possible to construct a literature model over a particular base.
Currently, this feature is only supported for toric base spaces.
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t2 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> length(singular_loci(t2))
2
```
Of course, this is also possible for Weierstrass models.
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, model_sections = Dict("b" => b), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> length(singular_loci(w))
1
```
For convenience, we also support a simplified constructor. Instead of the meta data of the article,
this constructor accepts an integer, which specifies the position of this model in our database.
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(3, base_space = B2, model_sections = Dict("b" => b), completeness_check = false)
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> length(singular_loci(w))
1
```
Similarly, also hypersurface models are supported:
```jldoctest
julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5")
Assuming that the first row of the given grading is the grading under Kbar

Hypersurface model over a not fully specified base

julia> explicit_model_sections(h)
Dict{String, MPolyRingElem} with 5 entries:
  "c2" => c2
  "c1" => c1
  "c3" => c3
  "b"  => b
  "c0" => c0

julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> h2 = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, model_sections = Dict("b" => b))
Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!

Hypersurface model over a concrete base

julia> hypersurface_equation_parametrization(h2)
b*w*v^2 - c0*u^4 - c1*u^3*v - c2*u^2*v^2 - c3*u*v^3 + w^2
```
"""
function literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
  model_dict = _find_model(doi, arxiv_id, version, equation)
  return literature_model(model_dict; model_parameters = model_parameters, base_space = base_space, model_sections = model_sections, completeness_check = completeness_check)
end

function literature_model(k::Int; model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
  model_dict = _find_model(k)
  return literature_model(model_dict; model_parameters = model_parameters, base_space = base_space, model_sections = model_sections, completeness_check = completeness_check)
end

function literature_model(model_dict::Dict{String, Any}; model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
  #return model_dict
  # (1) Deal with model parameters
  if haskey(model_dict, "model_parameters")
    needed_model_parameters = string.(model_dict["model_parameters"])
    
    # Make sure the user has provided values for all the model parameters
    for param in needed_model_parameters
      @req (param in keys(model_parameters)) "Some model parameters not provided; the given model requires these parameters:\n  $(join(needed_model_parameters, "\n  "))"
    end
    
    # Function to map a function of strings over arbitrarily nested Vectors of strings
    nested_string_map(f, s::String) = f(s)
    nested_string_map(f, v::Vector) = map(x -> nested_string_map(f, x), v)
    nested_string_map(f, a::Any) = a
    
    for (key, val) in model_parameters
      map!(x -> nested_string_map(s -> replace(s, key => string(val)), x), values(model_dict["model_data"]))
    end
  end
  
  
  # (2a) Construct the model over concrete base
  if dim(base_space) > 0
    
    # Currently, support only for toric bases
    @req base_space isa NormalToricVariety "Construction of literature models over concrete bases currently limited to toric bases"
    for (key, value) in model_sections
      @req value isa ToricDivisor "Construction of literature models over concrete bases currently requires toric divisors as model sections"
    end
    
    # Are all model sections specified?
    @req all(k->haskey(model_sections, k), model_dict["model_data"]["model_sections"]) "Not all model sections are specified"
    
    # Is the model specific for a base dimension? If so, make consistency check
    if haskey(model_dict["model_data"], "base_dim")
      @req dim(base_space) == Int(model_dict["model_data"]["base_dim"]) "Model requires base dimension different from dimension of provided base"
    end
    
    # Construct the model
    model = _construct_literature_model_over_concrete_base(model_dict, base_space, model_sections, completeness_check)
    @vprint :FTheoryModelPrinter 0 "Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!\n\n"
    
  # (2b) Construct the model over generic base
  else
    model = _construct_literature_model_over_arbitrary_base(model_dict)
  end
  
  
  # (3) Return the model after we set all required attributes
  _set_all_attributes(model, model_dict, model_parameters)
  return model
end



#######################################################
# 2. Helper function to find the specified model
#######################################################

function _find_model(doi::String, arxiv_id::String, version::String, equation::String)
  @req any(s -> s != "", [doi, arxiv_id, version, equation]) "No information provided; cannot perform look-up"
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  candidate_files = Vector{String}()
  for k in 1:length(file_index)
    if all([doi == "" || get(file_index[k], "journal_doi", nothing) == doi,
        arxiv_id == "" || get(file_index[k], "arxiv_id", nothing) == arxiv_id,
        version == "" || get(file_index[k], "arxiv_version", nothing) == version,
        equation == "" || get(file_index[k], "arxiv_equation", nothing) == equation])
      push!(candidate_files, string(file_index[k]["file"]))
    end
  end
  return _process_candidates(candidate_files)
end

function _find_model(l::Int)
  @req l >= 1 "Model index must be at least 1"
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  candidate_files = Vector{String}()
  for k in 1:length(file_index)
    if get(file_index[k], "model_index", nothing) == string(l)
      push!(candidate_files, string(file_index[k]["file"]))
    end
  end
  return _process_candidates(candidate_files)
end

function _process_candidates(candidate_files::Vector{String})
  @req length(candidate_files) != 0 "We could not find any models matching the given model index"
  @req(length(candidate_files) == 1,
    begin
      dicts = map(f -> JSON.parsefile(joinpath(@__DIR__, "Models/" * f)), candidate_files)
      dois = map(d -> get(d["journal_data"], "doi", nothing), dicts)
      ids = map(d -> get(d["arxiv_data"], "id", nothing), dicts)
      versions = map(d -> get(d["arxiv_data"], "version", nothing), dicts)
      equations = map(d -> get(d["arxiv_data"]["model_location"], "equation", nothing), dicts)
      strings = ["doi: $(dois[i]), arxiv_id: $(ids[i]), version: $(versions[i]), equation: $(equations[i])" for i in 1:length(dicts)]
      "We could not uniquely identify the model. The matched models have the following data:\n$(reduce((s1, s2) -> s1 * "\n" * s2, strings))"
    end)
  model_dict = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))
  model_dict["literature_identifier"] = candidate_files[1][6:end - 5]
  return model_dict
end


#######################################################
# 3. Constructing models over concrete bases
#######################################################

# Constructs literature model over concrete base
function _construct_literature_model_over_concrete_base(model_dict::Dict{String,Any}, base_space::FTheorySpace, model_sections::Dict{String,ToricDivisor}, completeness_check::Bool)

  # We first create a polynomial ring in which we can read all model sections as polynomials of the defining sections
  @req ((model_dict["model_descriptors"]["type"] == "tate") || (model_dict["model_descriptors"]["type"] == "weierstrass") || (model_dict["model_descriptors"]["type"] == "hypersurface")) "Model is not a Tate or Weierstrass model"
  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = polynomial_ring(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)
  vars = [string(g) for g in gens(auxiliary_base_ring)]

  # Make list of divisor classes which express the internal model sections.
  model_sections_divisor_list = vcat([anticanonical_divisor(base_space)], [model_sections[vars[k]] for k in 1:length(model_sections)])

  # Find divisor classes of the internal model sections
  auxiliary_base_grading = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in model_dict["model_data"]["auxiliary_base_grading"]]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i:i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  internal_model_sections = Dict{String, ToricDivisor}()
  for k in 1+length(model_sections):ngens(auxiliary_base_ring)
    divisor = sum([auxiliary_base_grading[l,k] * model_sections_divisor_list[l] for l in 1:nrows(auxiliary_base_grading)])
    internal_model_sections[vars[k]] = divisor
  end

  # Next, generate random values for all involved sections.
  explicit_model_sections = Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
  for (key, value) in model_sections
    @req is_effective(toric_divisor_class(value)) "Encountered a non-effective model section"
    #explicit_model_sections[key] = generic_section(toric_line_bundle(value));
    # Lead to error when computing singular loci - currently only monomials allowed...
    explicit_model_sections[key] = basis_of_global_sections(toric_line_bundle(value))[end]
  end
  for (key, value) in internal_model_sections
    @req is_effective(toric_divisor_class(value)) "Encountered a non-effective (internal) model section"
    explicit_model_sections[key] = generic_section(toric_line_bundle(value));
  end

  # Construct the model
  map = hom(auxiliary_base_ring, cox_ring(base_space), [explicit_model_sections[k] for k in vars])
  if model_dict["model_descriptors"]["type"] == "tate"

    # Compute Tate sections
    a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_base_ring)
    a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_base_ring)
    a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_base_ring)
    a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_base_ring)
    a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_base_ring)

    # Complete explicit_model_sections
    explicit_model_sections["a1"] = map(a1)
    explicit_model_sections["a2"] = map(a2)
    explicit_model_sections["a3"] = map(a3)
    explicit_model_sections["a4"] = map(a4)
    explicit_model_sections["a6"] = map(a6)

    # Find defining_section_parametrization
    defining_section_parametrization = Dict{String, MPolyRingElem}()
    if !("a1" in vars) || (a1 != eval_poly("a1", parent(a1)))
      defining_section_parametrization["a1"] = a1
    end
    if !("a2" in vars) || (a2 != eval_poly("a2", parent(a2)))
      defining_section_parametrization["a2"] = a2
    end
    if !("a3" in vars) || (a3 != eval_poly("a3", parent(a3)))
      defining_section_parametrization["a3"] = a3
    end
    if !("a4" in vars) || (a4 != eval_poly("a4", parent(a4)))
      defining_section_parametrization["a4"] = a4
    end
    if !("a6" in vars) || (a6 != eval_poly("a6", parent(a6)))
      defining_section_parametrization["a6"] = a6
    end

    # Create the model
    model = global_tate_model(base_space, explicit_model_sections, defining_section_parametrization; completeness_check = completeness_check)

  elseif model_dict["model_descriptors"]["type"] == "weierstrass"

    # Compute Weierstrass sections
    f = eval_poly(get(model_dict["model_data"], "f", "0"), auxiliary_base_ring)
    g = eval_poly(get(model_dict["model_data"], "g", "0"), auxiliary_base_ring)

    # Complete explicit_model_sections
    explicit_model_sections["f"] = map(f)
    explicit_model_sections["g"] = map(g)

    # Find defining_section_parametrization
    defining_section_parametrization = Dict{String, MPolyRingElem}()
    if !("f" in vars) || (f != eval_poly("f", parent(f)))
      defining_section_parametrization["f"] = f
    end
    if !("g" in vars) || (g != eval_poly("g", parent(g)))
      defining_section_parametrization["g"] = g
    end

    # Create the model
    model = weierstrass_model(base_space, explicit_model_sections, defining_section_parametrization; completeness_check = completeness_check)
  
  elseif model_dict["model_descriptors"]["type"] == "hypersurface"

    # Extract fiber ambient space
    rays = [[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_rays"]]
    max_cones = IncidenceMatrix([[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_max_cones"]])
    fas = normal_toric_variety(max_cones, rays; non_redundant = true)
    fiber_amb_coordinates = string.(model_dict["model_data"]["fiber_ambient_space_coordinates"])
    set_coordinate_names(fas, fiber_amb_coordinates)

    # Extract the divisor classes of the first two classes of the fiber...
    D1 = [a for a in model_dict["model_data"]["D1"]]
    D1_dc = toric_divisor_class(sum([D1[l] * model_sections_divisor_list[l] for l in 1:nrows(auxiliary_base_grading)]))
    D2 = [a for a in model_dict["model_data"]["D2"]]
    D2_dc = toric_divisor_class(sum([D2[l] * model_sections_divisor_list[l] for l in 1:nrows(auxiliary_base_grading)]))

    # Create the model
    model = hypersurface_model(base_space, fas, D1_dc, D2_dc; completeness_check = completeness_check)

    # Remember explicit model sections
    model.explicit_model_sections = explicit_model_sections

    # Remember hypersurface_parametrization
    auxiliary_ambient_ring, _ = polynomial_ring(QQ, vcat(vars, fiber_amb_coordinates), cached=false)
    parametrized_hypersurface_equation = eval_poly(model_dict["model_data"]["hypersurface_equation"], auxiliary_ambient_ring)
    model.hypersurface_equation_parametrization = parametrized_hypersurface_equation

    # Set explicit hypersurface equation
    images1 = [eval_poly(string(explicit_model_sections[k]), cox_ring(ambient_space(model))) for k in vars]
    images2 = [eval_poly(string(k), cox_ring(ambient_space(model))) for k in fiber_amb_coordinates]
    map = hom(parent(parametrized_hypersurface_equation), cox_ring(ambient_space(model)), vcat(images1, images2))
    model.hypersurface_equation = map(parametrized_hypersurface_equation)

  else

    @req false "Model is not a Tate, Weierstrass or hypersurface model"

  end

  # Return the model
  set_attribute!(model, :explicit_model_sections => explicit_model_sections)
  return model
end



#######################################################
# 4. Constructing models over arbitrary bases
#######################################################

# Constructs literature model over arbitrary base
function _construct_literature_model_over_arbitrary_base(model_dict::Dict{String,Any})
  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = polynomial_ring(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)

  @req haskey(model_dict["model_data"], "auxiliary_base_grading") "Database does not specify auxiliary_base_grading, but is vital for model constrution, so cannot proceed"
  auxiliary_base_grading = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in model_dict["model_data"]["auxiliary_base_grading"]]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i:i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  
  base_dim = get(model_dict["model_data"], "base_dim", 3)

  # Construct the model
  if model_dict["model_descriptors"]["type"] == "tate"

    a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_base_ring)
    a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_base_ring)
    a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_base_ring)
    a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_base_ring)
    a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_base_ring)
    model = global_tate_model(auxiliary_base_ring, auxiliary_base_grading, base_dim, [a1, a2, a3, a4, a6])

  elseif model_dict["model_descriptors"]["type"] == "weierstrass"

    f = eval_poly(get(model_dict["model_data"], "f", "0"), auxiliary_base_ring)
    g = eval_poly(get(model_dict["model_data"], "g", "0"), auxiliary_base_ring)
    model = weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, base_dim, f, g)

  elseif model_dict["model_descriptors"]["type"] == "hypersurface"

    # Extract base variable names
    auxiliary_base_vars = [string(g) for g in gens(auxiliary_base_ring)]

    # Extract fiber ambient space
    rays = [[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_rays"]]
    max_cones = IncidenceMatrix([[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_max_cones"]])
    fas = normal_toric_variety(max_cones, rays; non_redundant = true)
    fiber_amb_coordinates = string.(model_dict["model_data"]["fiber_ambient_space_coordinates"])
    set_coordinate_names(fas, fiber_amb_coordinates)

    # Extract the divisor classes of the first two classes of the fiber...
    D1 = [a for a in model_dict["model_data"]["D1"]]
    D2 = [a for a in model_dict["model_data"]["D2"]]

    # Extract the hypersurface equation
    ambient_ring, _ = polynomial_ring(QQ, vcat(auxiliary_base_vars, fiber_amb_coordinates), cached = false)
    p = eval_poly(model_dict["model_data"]["hypersurface_equation"], ambient_ring)
    
    # Create the model
    model = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, base_dim, fas, D1, D2, p)

  else

    @req false "Model is not a Tate, Weierstrass or hypersurface model"

  end
  return model
end



#######################################################
# 5. Functions for settings attributes of models
#######################################################

function _set_all_attributes(model::AbstractFTheoryModel, model_dict::Dict{String, Any}, model_parameters::Dict{String,<:Any})
  set_attribute!(model, :partially_resolved, false)
  set_literature_identifier(model, model_dict["literature_identifier"])
  
  set_model_description(model, model_dict["model_descriptors"]["description"])
  
  set_paper_authors(model, string.(model_dict["paper_metadata"]["authors"]))
  set_paper_buzzwords(model, string.(model_dict["paper_metadata"]["buzzwords"]))
  set_paper_description(model, model_dict["paper_metadata"]["description"])
  set_paper_title(model, model_dict["paper_metadata"]["title"])

  set_arxiv_doi(model, model_dict["arxiv_data"]["doi"])
  set_arxiv_link(model, model_dict["arxiv_data"]["link"])
  set_arxiv_id(model, model_dict["arxiv_data"]["id"])
  set_arxiv_version(model, model_dict["arxiv_data"]["version"])
  set_arxiv_model_equation_number(model, model_dict["arxiv_data"]["model_location"]["equation"])
  set_arxiv_model_page(model, model_dict["arxiv_data"]["model_location"]["page"])
  set_arxiv_model_section(model, model_dict["arxiv_data"]["model_location"]["section"])
  
  set_journal_doi(model, model_dict["journal_data"]["doi"])
  set_journal_link(model, model_dict["journal_data"]["link"])
  set_journal_year(model, model_dict["journal_data"]["year"])
  set_journal_volume(model, model_dict["journal_data"]["volume"])
  set_journal_name(model, model_dict["journal_data"]["journal"])
  if haskey(model_dict["journal_data"], "report_numbers")
    set_journal_report_numbers(model, string.(model_dict["journal_data"]["report_numbers"]))
  end
  set_journal_pages(model, model_dict["journal_data"]["pages"])
  set_journal_model_equation_number(model, model_dict["journal_data"]["model_location"]["equation"])
  set_journal_model_page(model, model_dict["journal_data"]["model_location"]["page"])
  set_journal_model_section(model, model_dict["journal_data"]["model_location"]["section"])
  
  if haskey(model_dict, "related_models")
    set_attribute!(model, :related_literature_models => [str[6:end - 5] for str in model_dict["related_models"]])
  end
  
  if haskey(model_dict, "associated_models")
    set_attribute!(model, :associated_literature_models => [str[6:end - 5] for str in model_dict["associated_models"]])
  end
  
  if haskey(model_dict, "model_parameters")
    set_attribute!(model, :model_parameters => model_parameters)
  end
  
  if haskey(model_dict["model_data"], "resolutions")
    set_resolutions(model, [[[string.(c) for c in r[1]], string.(r[2])] for r in model_dict["model_data"]["resolutions"]])
  end
  
  if haskey(model_dict["model_data"], "resolution_generating_sections")
    value = [[[string.(k) for k in sec] for sec in res] for res in model_dict["model_data"]["resolution_generating_sections"]]
    set_resolution_generating_sections(model, value)
  end
  
  if haskey(model_dict["model_data"], "resolution_zero_sections")
    set_resolution_zero_sections(model, [[string.(a) for a in b] for b in model_dict["model_data"]["resolution_zero_sections"]])
  end
  
  if haskey(model_dict["model_data"], "weighted_resolutions")
    set_weighted_resolutions(model, [[[[string.(c[1]), c[2]] for c in r[1]], string.(r[2])] for r in model_dict["model_data"]["weighted_resolutions"]])
  end
  
  if haskey(model_dict["model_data"], "weighted_resolution_generating_sections")
    value = [[[string.(k) for k in sec] for sec in res] for res in model_dict["model_data"]["weighted_resolution_generating_sections"]]
    set_weighted_resolution_generating_sections(model, value)
  end
  
  if haskey(model_dict["model_data"], "weighted_resolution_zero_sections")
    set_weighted_resolution_zero_sections(model, [[string.(a) for a in b] for b in model_dict["model_data"]["weighted_resolution_zero_sections"]])
  end
  
  if haskey(model_dict["model_data"], "zero_section")
    set_zero_section(model, string.(model_dict["model_data"]["zero_section"]))
  end

  if haskey(model_dict["model_data"], "generating_sections")
    set_generating_sections(model, map(k -> string.(k), model_dict["model_data"]["generating_sections"]))
  end
end



#######################################################
# 6. Function to display all known literature models
#######################################################

function display_all_literature_models()
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  sorted_dicts = sort(file_index, by = x -> parse(Int, x["model_index"]))
  for dict in sorted_dicts
    print("Model $(dict["model_index"]):\n")
    print(dict)
    print("\n\n")
  end
end
