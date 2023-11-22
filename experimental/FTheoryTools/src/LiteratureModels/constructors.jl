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
Normal toric variety

julia> a1,a21,a32,a43,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = domain(blow_up(v, I))
Normal toric variety

julia> cox_ring(v2)
Multivariate polynomial ring in 9 variables over QQ graded by
  a1 -> [1 0 0 0]
  a21 -> [0 1 0 0]
  a32 -> [-1 2 0 0]
  a43 -> [-2 3 0 0]
  w -> [0 0 1 0]
  x -> [0 1 1 2]
  y -> [1 1 1 3]
  z -> [0 0 0 1]
  e -> [2 -1 -1 0]
```
It is also possible to construct a literature model over a particular base.
Currently, this feature is only supported for toric base spaces.
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t2 = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, model_sections = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> length(singular_loci(t2))
2
```
"""
function literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
  
  # (1) Find the model
  model_dict = _find_model(doi, arxiv_id, version, equation)
  
  
  # (2) Deal with model parameters
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
  
  
  # (3a) Construct the model over concrete base
  if dim(base_space) > 0
    
    # Currently, support only for toric bases
    @req typeof(base_space) == NormalToricVariety "Construction of literature models over concrete bases currently limited to toric bases"
    for (key, value) in model_sections
      @req typeof(value) == ToricDivisor "Construction of literature models over concrete bases currently requires toric divisors as model sections"
    end
    
    # Are all model sections specified?
    needed_model_sections = model_dict["model_data"]["model_sections"]
    @req all(k->haskey(model_sections, k), needed_model_sections) "Not all model sections are specified"
    
    # Is the model specific for a base dimension? If so, make consistency check
    if haskey(model_dict["model_data"], "base_dim")
      @req dim(base_space) == Int(model_dict["model_data"]["base_dim"]) "Model requires base dimension different from dimension of provided base"
    end
    
    # Appropriately construct the model
    if model_dict["model_descriptors"]["type"] == "tate"
      model = _construct_literature_tate_model(model_dict, base_space, model_sections, completeness_check)
    elseif model_dict["model_descriptors"]["type"] == "weierstrass"
      @req false "Weierstrass literature models can currently not be created over a concrete base space"
    else
      @req false "Model is not a Tate or Weierstrass model"
    end
    
  # (3b) Construct the model over a generic base
  else
    if model_dict["model_descriptors"]["type"] == "tate"
      model = _construct_literature_tate_model(model_dict)
    elseif model_dict["model_descriptors"]["type"] == "weierstrass"
      model = _construct_literature_weierstrass_model(model_dict)
    else
      @req false "Model is not a Tate or Weierstrass model"
    end
  end

  # Return the model
  _set_all_attributes(model, model_dict, model_parameters)
  return model
end


#######################################################
# 2. Helper function to find the specified model
#######################################################

function _find_model(doi::String, arxiv_id::String, version::String, equation::String)

  # Check that we have at least some information...
  @req any(s -> s != "", [doi, arxiv_id, version, equation]) "No information provided; cannot perform look-up"

  # Create list of possible candidate files
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

  # Check if we found exactly one file, i.e., we were able to identify the model uniquely
  @req length(candidate_files) != 0 "We could not find any models matching the given identifiers"
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

  # Create dictionary
  model_dict = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))
  # Add literature identifier. For this, remove 'model' from the front and '.json' from the end of the file name.
  model_dict["literature_identifier"] = candidate_files[1][6:end - 5]

  # Return the dictionary
  return model_dict
end



#######################################################
# 3. Constructing models over concrete bases
#######################################################

# Constructs Tate model over concrete base from given Tate literature model
function _construct_literature_tate_model(model_dict::Dict{String,Any}, base_space::FTheorySpace, model_sections::Dict{String,ToricDivisor}, completeness_check::Bool)

  # Generate random sections of the model sections
  explicit_model_sections = Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
  for (key, value) in model_sections
    # Lead to an error when trying to compute the singular loci.
    # Likely, the current implementation assumes that the singular locus is given by the vanishing of a single monom?
    explicit_model_sections[key] = basis_of_global_sections(toric_line_bundle(value))[end]
    #explicit_model_sections[key] = sum([rand(Int) * b for b in basis_of_global_sections(toric_line_bundle(value))]);
  end

  # For a Tate model we need to create a1, a2, a3, a4, a6 with ai a section of Kbar^i.
  # These are parametrized by the model sections above.
  # We need to extract the parametrization from the model_dict and then identify the additional parameters' toric divisors.

  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = polynomial_ring(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)
  a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_base_ring)
  a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_base_ring)
  a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_base_ring)
  a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_base_ring)
  a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_base_ring)

  # Currently, we assume that Tate sections are monomial expressions in the model sections and remaining parameters
  @req all(k->length(monomials(k)) == 1 || is_zero(k), [a1, a2, a3, a4, a6]) "Tate sections are assumed to be monomials in model sections and remaining parameters"

  # Initialize necessary structures
  divisor_list = Vector{ToricDivisor}()
  coefficient_matrix = Vector{Vector{Int64}}()
  Kbar = anticanonical_divisor(base_space)
  cut_off = ngens(auxiliary_base_ring) - length(explicit_model_sections)
  variables = gens(auxiliary_base_ring)
  model_section_divisor_vector = Vector{ToricDivisor}()
  for k in cut_off+1:length(variables)
    push!(model_section_divisor_vector, model_sections[string(variables[k])])
  end

  # Fill divisor_list and coefficient_matrix
  if is_zero(a1) == false
    exp_vector = exponent_vector(a1,1)
    model_section_combination = sum([exp_vector[k] * model_section_divisor_vector[k - cut_off] for k in cut_off+1:length(exp_vector)])
    push!(divisor_list, Kbar - model_section_combination)
    push!(coefficient_matrix, [exp_vector[k] for k in 1:cut_off])
  end
  if is_zero(a2) == false
    exp_vector = exponent_vector(a2,1)
    model_section_combination = sum([exp_vector[k] * model_section_divisor_vector[k - cut_off] for k in cut_off+1:length(exp_vector)])
    push!(divisor_list, 2 * Kbar - model_section_combination)
    push!(coefficient_matrix, [exp_vector[k] for k in 1:cut_off])
  end
  if is_zero(a3) == false
    exp_vector = exponent_vector(a3,1)
    model_section_combination = sum([exp_vector[k] * model_section_divisor_vector[k - cut_off] for k in cut_off+1:length(exp_vector)])
    push!(divisor_list, 3 * Kbar - model_section_combination)
    push!(coefficient_matrix, [exp_vector[k] for k in 1:cut_off])
  end
  if is_zero(a4) == false
    exp_vector = exponent_vector(a4,1)
    model_section_combination = sum([exp_vector[k] * model_section_divisor_vector[k - cut_off] for k in cut_off+1:length(exp_vector)])
    push!(divisor_list, 4 * Kbar - model_section_combination)
    push!(coefficient_matrix, [exp_vector[k] for k in 1:cut_off])
  end
  if is_zero(a6) == false
    exp_vector = exponent_vector(a6,1)
    model_section_combination = sum([exp_vector[k] * model_section_divisor_vector[k - cut_off] for k in cut_off+1:length(exp_vector)])
    push!(divisor_list, 6 * Kbar - model_section_combination)
    push!(coefficient_matrix, [exp_vector[k] for k in 1:cut_off])
  end

  # Find a left inverse of the coefficient matrix
  inverse_matrix = solve_left(matrix(ZZ, coefficient_matrix), identity_matrix(ZZ, cut_off))

  # Now we are in the position to tell the toric divisors for the internal parameters
  internal_model_sections = Dict{String, ToricDivisor}()
  for k in 1:cut_off
    coeffs = inverse_matrix[k,:]
    divisor = sum([coeffs[l] * divisor_list[l] for l in 1:length(coeffs)])
    internal_model_sections[string(variables[k])] = divisor
  end

  # Compute random internal model sections
  for (key, value) in internal_model_sections
    explicit_model_sections[key] = sum([rand(Int) * b for b in basis_of_global_sections(toric_line_bundle(value))]);
  end

  # Finally, let us create the actual Tate sections
  explicit_a1 = zero(cox_ring(base_space))
  if is_zero(a1) == false
    exp_vector = exponent_vector(a1,1)
    explicit_a1 = prod([explicit_model_sections[string(variables[k])]^exp_vector[k] for k in 1:ngens(auxiliary_base_ring)])
  end
  explicit_a2 = zero(cox_ring(base_space))
  if is_zero(a2) == false
    exp_vector = exponent_vector(a2,1)
    explicit_a2 = prod([explicit_model_sections[string(variables[k])]^exp_vector[k] for k in 1:ngens(auxiliary_base_ring)])
  end
  explicit_a3 = zero(cox_ring(base_space))
  if is_zero(a3) == false
    exp_vector = exponent_vector(a3,1)
    explicit_a3 = prod([explicit_model_sections[string(variables[k])]^exp_vector[k] for k in 1:ngens(auxiliary_base_ring)])
  end
  explicit_a4 = zero(cox_ring(base_space))
  if is_zero(a4) == false
    exp_vector = exponent_vector(a4,1)
    explicit_a4 = prod([explicit_model_sections[string(variables[k])]^exp_vector[k] for k in 1:ngens(auxiliary_base_ring)])
  end
  explicit_a6 = zero(cox_ring(base_space))
  if is_zero(a6) == false
    exp_vector = exponent_vector(a6,1)
    explicit_a6 = prod([explicit_model_sections[string(variables[k])]^exp_vector[k] for k in 1:ngens(auxiliary_base_ring)])
  end

  # Construct the model
  ais = [explicit_a1, explicit_a2, explicit_a3, explicit_a4, explicit_a6]
  model = global_tate_model(base_space, ais; completeness_check = completeness_check)

  # Remember the explicit model sections for autoresolution
  set_attribute!(model, :explicit_model_sections => explicit_model_sections)

  # Return the model
  return model
end



#######################################################
# 4. Constructing models over arbitrary bases
#######################################################

# Constructs Tate model from given Tate literature model
function _construct_literature_tate_model(model_dict::Dict{String,Any})
  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = polynomial_ring(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)
  
  base_dim = get(model_dict["model_data"], "base_dim", 3)
  
  a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_base_ring)
  a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_base_ring)
  a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_base_ring)
  a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_base_ring)
  a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_base_ring)
  
  @req haskey(model_dict["model_data"], "auxiliary_base_grading") "Database does not specify auxiliary_base_grading, but is vital for model constrution, so cannot proceed"
  auxiliary_base_grading = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in model_dict["model_data"]["auxiliary_base_grading"]]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  return global_tate_model(auxiliary_base_ring, auxiliary_base_grading, base_dim, [a1, a2, a3, a4, a6])
end

# Constructs Weierstrass model from given Weierstrass literature model
function _construct_literature_weierstrass_model(model_dict::Dict{String,Any})
  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = polynomial_ring(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)
  
  base_dim = get(model_dict["model_data"], "base_dim", 3)
  
  f = eval_poly(get(model_dict["model_data"], "f", "0"), auxiliary_base_ring)
  g = eval_poly(get(model_dict["model_data"], "g", "0"), auxiliary_base_ring)
  
  @req haskey(model_dict["model_data"], "auxiliary_base_grading") "Database does not specify auxiliary_base_grading, but is vital for model constrution, so cannot proceed"
  auxiliary_base_grading = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in model_dict["model_data"]["auxiliary_base_grading"]]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  return weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, base_dim, f, g)
end




#######################################################
# 5. Functions for settings attributes of models
#######################################################

function _set_all_attributes(model::AbstractFTheoryModel, model_dict::Dict{String, Any}, model_parameters::Dict{String,<:Any})
  set_attribute!(model, :literature_identifier => model_dict["literature_identifier"])
  set_attribute!(model, :partially_resolved, false)
  
  _set_model_attribute(model, model_dict, "model_descriptors", "description", "model_description")
  
  _set_model_attribute(model, model_dict, "paper_metadata", "authors", "paper_authors")
  _set_model_attribute(model, model_dict, "paper_metadata", "buzzwords", "paper_buzzwords")
  _set_model_attribute(model, model_dict, "paper_metadata", "description", "paper_description")
  _set_model_attribute(model, model_dict, "paper_metadata", "title", "paper_title")
  
  _set_model_attribute(model, model_dict, "arxiv_data", "doi", "arxiv_doi")
  _set_model_attribute(model, model_dict, "arxiv_data", "link", "arxiv_link")
  _set_model_attribute(model, model_dict, "arxiv_data", "id", "arxiv_id")
  _set_model_attribute(model, model_dict, "arxiv_data", "version", "arxiv_version")  
  _set_model_attribute(model, model_dict, "arxiv_data", "model_location", "equation", "arxiv_model_equation_number")
  _set_model_attribute(model, model_dict, "arxiv_data", "model_location", "page", "arxiv_model_page")
  _set_model_attribute(model, model_dict, "arxiv_data", "model_location", "section", "arxiv_model_section")
  
  _set_model_attribute(model, model_dict, "journal_data", "doi", "journal_doi")
  _set_model_attribute(model, model_dict, "journal_data", "link", "journal_link")
  _set_model_attribute(model, model_dict, "journal_data", "year", "journal_year")
  _set_model_attribute(model, model_dict, "journal_data", "volume", "journal_volume")
  _set_model_attribute(model, model_dict, "journal_data", "report_numbers", "journal_report_numbers")
  _set_model_attribute(model, model_dict, "journal_data", "pages", "journal_pages")
  _set_model_attribute(model, model_dict, "journal_data", "model_location", "equation", "journal_model_equation_number")
  _set_model_attribute(model, model_dict, "journal_data", "model_location", "page", "journal_model_page")
  _set_model_attribute(model, model_dict, "journal_data", "model_location", "section", "journal_model_section")
  
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
    resolutions_data = model_dict["model_data"]["resolutions"]
    resolutions = [[[string.(center) for center in res[1]], string.(res[2])] for res in resolutions_data]
    set_attribute!(model, :resolutions => resolutions)
  end
  
  if haskey(model_dict["model_data"], "resolution_generating_sections")
    resolution_generating_sections_data = model_dict["model_data"]["resolution_generating_sections"]
    resolution_generating_sections = [[[string.(factor) for factor in sec] for sec in res] for res in resolution_generating_sections_data]
    set_attribute!(model, :resolution_generating_sections => resolution_generating_sections)
  end
  
  if haskey(model_dict["model_data"], "resolution_zero_sections")
    resolution_zero_sections_data = model_dict["model_data"]["resolution_zero_sections"]
    resolution_zero_sections = [[string.(factor) for factor in res] for res in resolution_zero_sections_data]
    set_attribute!(model, :resolution_zero_sections => resolution_zero_sections)
  end
  
  if haskey(model_dict["model_data"], "weighted_resolutions")
    weighted_resolutions_data = model_dict["model_data"]["weighted_resolutions"]
    weighted_resolutions = [[[[string.(center[1]), center[2]] for center in res[1]], string.(res[2])] for res in weighted_resolutions_data]
    set_attribute!(model, :weighted_resolutions => weighted_resolutions)
  end
  
  if haskey(model_dict["model_data"], "weighted_resolution_generating_sections")
    weighted_resolution_generating_sections_data = model_dict["model_data"]["weighted_resolution_generating_sections"]
    weighted_resolution_generating_sections = [[[string.(factor) for factor in sec] for sec in res] for res in weighted_resolution_generating_sections_data]
    set_attribute!(model, :weighted_resolution_generating_sections => weighted_resolution_generating_sections)
  end
  
  if haskey(model_dict["model_data"], "weighted_resolution_zero_sections")
    weighted_resolution_zero_sections_data = model_dict["model_data"]["weighted_resolution_zero_sections"]
    weighted_resolution_zero_sections = [[string.(factor) for factor in res] for res in weighted_resolution_zero_sections_data]
    set_attribute!(model, :weighted_resolution_zero_sections => weighted_resolution_zero_sections)
  end
  
  if typeof(model.base_space) == NormalToricVariety
    base_ring = cox_ring(model.base_space) # THIS CURRENTLY ASSUMES THE BASE IS TORIC, SHOULD FIX
    if haskey(model_dict["model_data"], "zero_section")
      set_attribute!(model, :zero_section => [eval_poly(coord, base_ring) for coord in model_dict["model_data"]["zero_section"]])
    end
    if haskey(model_dict["model_data"], "generating_sections")
      set_attribute!(model, :generating_sections => [[eval_poly(coord, base_ring) for coord in gen_sec] for gen_sec in model_dict["model_data"]["generating_sections"]])
    end
  end
  
end

function _set_model_attribute(m::AbstractFTheoryModel, m_dict::Dict{String, Any}, l::String, t::String, t_name::String)
  if haskey(m_dict[l], t)
    if typeof(m_dict[l][t]) <: Vector{Any}
      set_attribute!(m, Symbol(t_name) => [String(k) for k in m_dict[l][t]])
    else
      set_attribute!(m, Symbol(t_name) => m_dict[l][t])
    end
  end
end

function _set_model_attribute(m::AbstractFTheoryModel, m_dict::Dict{String, Any}, l1::String, l2::String, t::String, t_name::String)
  if haskey(m_dict[l1][l2], t)
    if typeof(m_dict[l1][l2][t]) <: Vector{Any}
      set_attribute!(m, Symbol(t_name) => [String(k) for k in m_dict[l1][l2][t]])
    else
      set_attribute!(m, Symbol(t_name) => m_dict[l1][l2][t])
    end
  end
end
