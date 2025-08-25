#######################################################
# 1. User interface for literature models
#######################################################

@doc raw"""
    literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), defining_classes::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)

Return a model from the F-theory literature identified by bibliographic metadata such as arXiv ID, DOI, or equation reference.
Many such models are well-known in the field—e.g., the *U(1)-restricted SU(5)-GUT model
[Krause, Mayrhofer, Weigand 2011](@cite KMW12)*—and are stored in a curated database accessible through this method.

To identify and retrieve a model from the database, supply one or more of the following optional arguments:

* `doi`: A string representing the DOI of the publication in which the model was introduced.
* `arxiv_id`: A string specifying the arXiv identifier of the preprint version or the journal print version.
* `version`: A string specifying the arXiv version (e.g. `"v1"`).
* `equation`: A string indicating the equation label used to define the model within the source publication.

The method attempts to match the given identifiers with a unique entry in the database. If no match is found, or if the information
is ambiguous (i.e., multiple matches exist), an error is raised.

Some literature models require additional input to be uniquely determined and constructed. For such cases, more optional arguments must be provided:

* `model_parameters`: A dictionary specifying parameters needed to fix a model within a family, e.g. `Dict("k" => 5)`.
* `defining_classes`: A dictionary specifying divisor classes necessary to fully define the geometry.
* `base_space`: Optionally specify a concrete base over which the model should be constructed. Currently, only toric base spaces are supported.
* `completeness_check`: Set this to `false` to skip time consuming completeness checks of the base geometry to gain more performance.

See the [Literature Models](@ref literature_models) documentation page for a deeper discussion of these fields.

First, notice how you can create the global Tate model from [Krause, Mayrhofer, Weigand 2011](@cite KMW12)
over an unspecified base:

# Examples
```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> coordinate_ring(ambient_space(t))
Multivariate polynomial ring in 8 variables over QQ graded by
  w -> [0 1 0]
  a1 -> [1 0 0]
  a21 -> [2 -1 0]
  a32 -> [3 -2 0]
  a43 -> [4 -3 0]
  x -> [2 0 2]
  y -> [3 0 3]
  z -> [0 0 1]
```

Certainly, this is also possible over a concrete base, provided that the defining classes are provided:

# Examples
```jldoctest
julia> B3 = projective_space(NormalToricVariety, 3)
Normal toric variety

julia> w = torusinvariant_prime_divisors(B3)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1", base_space = B3, defining_classes = Dict("w" => w), completeness_check = false)
Global Tate model over a concrete base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)
```

Certainly, this can be mimiced for Weierstrass models, e.g. [Morrison, Park 2012](@cite MP12):

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(arxiv_id = "1208.2695", equation = "B.19", base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)
Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> using Random;

julia> length(singular_loci(w; rng = Random.Xoshiro(1234)))
1
```

For convenience, we also support a simplified constructor. Instead of the meta data of the article,
this constructor accepts an integer, which specifies the position of this model in our database.
Here is how this can be used to create the model [Morrison, Park 2012](@cite MP12):

# Examples
```jldoctest
julia> B2 = projective_space(NormalToricVariety, 2)
Normal toric variety

julia> b = torusinvariant_prime_divisors(B2)[1]
Torus-invariant, prime divisor on a normal toric variety

julia> w = literature_model(3, base_space = B2, defining_classes = Dict("b" => b), completeness_check = false)
Weierstrass model over a concrete base -- U(1) Weierstrass model based on arXiv paper 1208.2695 Eq. (B.19)

julia> using Random;

julia> length(singular_loci(w; rng = Random.Xoshiro(1234)))
1
```

Finally, let us showcase that we also support hypersurface model constructions from the literature.
For this, we stay with [Morrison, Park 2012](@cite MP12), but this time create this model as a
hypersurface model. This works as follows:

# Examples
```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> h = literature_model(arxiv_id = "1208.2695", equation = "B.5")
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

julia> h2 = literature_model(arxiv_id = "1208.2695", equation = "B.5", base_space = B2, defining_classes = Dict("b" => b))
Hypersurface model over a concrete base

julia> hypersurface_equation_parametrization(h2)
b*w*v^2 - c0*u^4 - c1*u^3*v - c2*u^2*v^2 - c3*u*v^3 + w^2
```
"""
function literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="", type::String="", model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), defining_classes::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
  model_dict = _find_model(doi, arxiv_id, version, equation, type)
  return literature_model(model_dict; model_parameters = model_parameters, base_space = base_space, model_sections = model_sections, defining_classes = defining_classes, completeness_check = completeness_check)
end

function literature_model(k::Int; model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String,Any}(), defining_classes::Dict{String, <:Any} = Dict{String,Any}(), completeness_check::Bool = true)
  model_dict = _find_model(k)
  return literature_model(model_dict; model_parameters = model_parameters, base_space = base_space, model_sections = model_sections, defining_classes = defining_classes, completeness_check = completeness_check)
end

function literature_model(model_dict::Dict{String, Any}; model_parameters::Dict{String,<:Any} = Dict{String,Any}(), base_space::FTheorySpace = affine_space(NormalToricVariety, 0), model_sections::Dict{String, <:Any} = Dict{String, Any}(), defining_classes::Dict{String, <:Any} = Dict{String, Any}(), completeness_check::Bool = true)
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
      map!(x -> nested_string_map(s -> replace(s, string("#", key) => string(val)), x), values(model_dict["model_descriptors"]))
      map!(x -> nested_string_map(s -> replace(s, r"\(([^(),]+)\)" => dim -> string("(", Oscar.eval_poly(string.(match(r"\(([^(),]+)\)", dim).captures[1]), ZZ),")")), x), values(model_dict["model_descriptors"]))
    end
  end
  
  # (2) The QSM need special treatment...
  if model_dict["arxiv_data"]["id"] == "1903.00009"
    model_dict["literature_identifier"] = "1903_00009"
    k = model_parameters["k"]
    qsmd_path = artifact"QSMDB"
    qsm_model = load(joinpath(qsmd_path, "$k.mrdi"))
    return qsm_model
  end

  # (2b) The F-theory model with the largest number of flux vacua needs special attention
  if model_dict["arxiv_data"]["id"] == "1511.03209"

    model_data_path = artifact"FTM-1511-03209/1511-03209.mrdi"
    return load(model_data_path)

    # Old code to create this model from scratch. I leave this here, so we can go back if needed.
    #=
    directory = joinpath(@__DIR__, "Models/1511_03209/1511-03209-base-space.mrdi")
    base_space = load(directory)
    set_attribute!(base_space, :coordinate_names, ["w$i" for i in 0:100])
    model = global_tate_model(base_space, completeness_check = false)
    model_dict["literature_identifier"] = "1511_03209"
    _set_all_attributes(model, model_dict, model_parameters)
    return model
    =#
    
  end

  # (3) Construct the model over concrete or arbitrary base
  if dim(base_space) > 0
    
    # Currently, support only for toric bases
    @req base_space isa NormalToricVariety "Construction of literature models over concrete bases currently limited to toric bases"

    # FIXME: Append model_sections to defining_classes. This might need fixing/extending in the future...
    defining_classes_provided = merge(model_sections, defining_classes)
    for (key, value) in defining_classes_provided
      if (value isa ToricDivisor || value isa ToricDivisorClass || value isa ToricLineBundle) == false
        error("Construction of literature models over concrete bases currently requires defining classes (and model sections) to be provided as toric divisor (classes) or line bundles")
      end
    end
    @req all(k->haskey(defining_classes_provided, k), model_dict["model_data"]["defining_classes"]) "Not all defining classes are specified"
    
    # Is the model specific for a base dimension? If so, make consistency check
    if haskey(model_dict["model_data"], "base_dim")
      @req dim(base_space) == Int(model_dict["model_data"]["base_dim"]) "Model requires base dimension different from dimension of provided base"
    end
    
    # Add additional information that is always known for weierstrass/tate 
    if (model_dict["model_descriptors"]["type"] == "weierstrass") || (model_dict["model_descriptors"]["type"] == "tate")
      model_dict["model_data"]["zero_section_class"] = "z"
    end

    # Construct the model
    model = _construct_literature_model_over_concrete_base(model_dict, base_space, defining_classes_provided, completeness_check)
    @vprint :FTheoryModelPrinter 1 "Construction over concrete base may lead to singularity enhancement. Consider computing singular_loci. However, this may take time!\n\n"
    
  else
    model = _construct_literature_model_over_arbitrary_base(model_dict)
  end
  
  
  # (4) Return the model after we set all required attributes
  _set_all_attributes(model, model_dict, model_parameters)
  return model
end



#######################################################
# 2. Helper function to find the specified model
#######################################################

function _find_model(doi::String, arxiv_id::String, version::String, equation::String, type::String)
  @req any(!isempty, [doi, arxiv_id, version, equation]) "No information provided; cannot perform look-up"
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  candidate_files = Vector{String}()
  for k in 1:length(file_index)
    if all([doi == "" || get(file_index[k], "journal_doi", nothing) == doi,
        arxiv_id == "" || get(file_index[k], "arxiv_id", nothing) == arxiv_id,
        version == "" || get(file_index[k], "arxiv_version", nothing) == version,
        equation == "" || get(file_index[k], "arxiv_equation", nothing) == equation,
        type == "" || get(file_index[k], "type", nothing) == type])
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
      types = map(d -> get(d["model_descriptors"], "type", nothing), dicts)
      strings = ["doi: $(dois[i]), arxiv_id: $(ids[i]), version: $(versions[i]), equation: $(equations[i]), type: $(types[i])" for i in 1:length(dicts)]
      "We could not uniquely identify the model. The matched models have the following data:\n$(reduce((s1, s2) -> s1 * "\n" * s2, strings))"
    end)
  model_dict = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))
  model_dict["literature_identifier"] = candidate_files[1][6:end - 5]
  return model_dict
end


#######################################################
# 3. Constructing models over concrete bases
#######################################################

# Construct literature model over concrete base
function _construct_literature_model_over_concrete_base(model_dict::Dict{String,Any}, base_space::FTheorySpace, defining_classes::Dict{String, <:Any}, completeness_check::Bool)

  # Make list and dict of the defining divisor classes
  defng_cls = string.(model_dict["model_data"]["defining_classes"])
  Kbar_and_defng_cls_as_divisor_classes = [anticanonical_divisor_class(base_space)]
  defining_classes_of_model = Dict()
  for k in 1:length(defng_cls)
    class_candidate = defining_classes[defng_cls[k]]
    if class_candidate isa ToricDivisor || class_candidate isa ToricLineBundle
      class_candidate = toric_divisor_class(class_candidate)
    end
    push!(Kbar_and_defng_cls_as_divisor_classes, class_candidate)
    defining_classes_of_model[defng_cls[k]] = class_candidate
  end

  # Find divisor classes of all tunable sections.
  @req haskey(model_dict["model_data"], "tunable_sections") "Database does not specify model sections for given model"
  tune_sec_names = string.(model_dict["model_data"]["tunable_sections"])
  cfs = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in model_dict["model_data"]["classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes"]]...)))
  cfs = vcat([[Int(k) for k in cfs[i:i,:]] for i in 1:nrows(cfs)]...)
  model_dict["model_data"]["classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes"] = cfs
  cl_of_secs = Dict(tune_sec_names[k] => sum(cfs[l, k] * Kbar_and_defng_cls_as_divisor_classes[l] for l in 1:nrows(cfs)) for k in 1:length(tune_sec_names))

  # Next, generate random values for all involved sections.
  model_sections = Dict{String, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}()
  for (key, value) in cl_of_secs
    @req is_effective(value) "Encountered a non-effective (internal) divisor class"
    if model_dict["arxiv_data"]["id"] == "1109.3454" && key == "w" && dim(base_space) == 3
      if torsion_free_rank(class_group_with_map(base_space)[1]) == 1 && degree(toric_line_bundle(cl_of_secs["w"])) == 1
        model_sections[key] = basis_of_global_sections(toric_line_bundle(value))[end]
      else
        model_sections[key] = generic_section(toric_line_bundle(value))
      end
    else
      model_sections[key] = generic_section(toric_line_bundle(value))
    end
  end

  # Construct the model
  auxiliary_ring, _ = polynomial_ring(QQ, tune_sec_names, cached=false)
  map = hom(auxiliary_ring, coordinate_ring(base_space), [model_sections[k] for k in tune_sec_names])
  if model_dict["model_descriptors"]["type"] == "tate"

    # Compute Tate sections
    a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_ring)
    a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_ring)
    a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_ring)
    a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_ring)
    a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_ring)

    # Compute defining model sections
    model_sections["a1"] = map(a1)
    model_sections["a2"] = map(a2)
    model_sections["a3"] = map(a3)
    model_sections["a4"] = map(a4)
    model_sections["a6"] = map(a6)

    # Find model_section_parametrization
    model_section_parametrization = Dict{String, MPolyRingElem}()
    if !("a1" in tune_sec_names) || (a1 != eval_poly("a1", parent(a1)))
      model_section_parametrization["a1"] = a1
    end
    if !("a2" in tune_sec_names) || (a2 != eval_poly("a2", parent(a2)))
      model_section_parametrization["a2"] = a2
    end
    if !("a3" in tune_sec_names) || (a3 != eval_poly("a3", parent(a3)))
      model_section_parametrization["a3"] = a3
    end
    if !("a4" in tune_sec_names) || (a4 != eval_poly("a4", parent(a4)))
      model_section_parametrization["a4"] = a4
    end
    if !("a6" in tune_sec_names) || (a6 != eval_poly("a6", parent(a6)))
      model_section_parametrization["a6"] = a6
    end

    # Create the model
    model = global_tate_model(base_space, model_sections, model_section_parametrization; completeness_check = completeness_check)

  elseif model_dict["model_descriptors"]["type"] == "weierstrass"

    # Compute Weierstrass sections
    f = eval_poly(get(model_dict["model_data"], "f", "0"), auxiliary_ring)
    g = eval_poly(get(model_dict["model_data"], "g", "0"), auxiliary_ring)

    # Compute defining model sections
    model_sections["f"] = map(f)
    model_sections["g"] = map(g)

    # Find model_section_parametrization
    model_section_parametrization = Dict{String, MPolyRingElem}()
    if !("f" in tune_sec_names) || (f != eval_poly("f", parent(f)))
      model_section_parametrization["f"] = f
    end
    if !("g" in tune_sec_names) || (g != eval_poly("g", parent(g)))
      model_section_parametrization["g"] = g
    end

    # Create the model
    model = weierstrass_model(base_space, model_sections, model_section_parametrization; completeness_check = completeness_check)

  elseif model_dict["model_descriptors"]["type"] == "hypersurface"

    # Extract fiber ambient space
    rays = [[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_rays"]]
    max_cones = IncidenceMatrix([[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_max_cones"]])
    fas = normal_toric_variety(max_cones, rays; non_redundant = true)
    fiber_amb_coordinates = string.(model_dict["model_data"]["fiber_ambient_space_coordinates"])
    set_coordinate_names(fas, fiber_amb_coordinates)

    # Extract the base divisor classes of the fiber coordinates
    fiber_twist_matrix = transpose(matrix(ZZ, (hcat(model_dict["model_data"]["fiber_twist_matrix"]...))))
    @req ncols(fiber_twist_matrix) == length(fiber_amb_coordinates) "Number of fiber coordinate names does not match number of provided fiber gradings"
    @req ncols(fiber_twist_matrix) == n_rays(fas) "Number of rays does not match number of provided fiber gradings"
    fiber_twist_divisor_classes = [sum([fiber_twist_matrix[l, k] * Kbar_and_defng_cls_as_divisor_classes[l] for l in 1:nrows(fiber_twist_matrix)]) for k in 1:ncols(fiber_twist_matrix)]

    # Compute the hypersurface equation and its parametrization
    auxiliary_ambient_ring, _ = polynomial_ring(QQ, vcat(tune_sec_names, fiber_amb_coordinates), cached=false)
    parametrized_hypersurface_equation = eval_poly(model_dict["model_data"]["hypersurface_equation"], auxiliary_ambient_ring)
    base_coordinates = string.(gens(coordinate_ring(base_space)))
    auxiliary_ambient_ring2, _ = polynomial_ring(QQ, vcat(base_coordinates, fiber_amb_coordinates), cached=false)
    images1 = [eval_poly(string(model_sections[k]), auxiliary_ambient_ring2) for k in tune_sec_names]
    images2 = [eval_poly(string(k), auxiliary_ambient_ring2) for k in fiber_amb_coordinates]
    map = hom(auxiliary_ambient_ring, auxiliary_ambient_ring2, vcat(images1, images2))
    hyper_equ = map(parametrized_hypersurface_equation)

    # Create the model
    model = hypersurface_model(base_space, fas, fiber_twist_divisor_classes, hyper_equ; completeness_check = completeness_check)
    model.hypersurface_equation_parametrization = parametrized_hypersurface_equation

  else
    error("Model is not a Tate, Weierstrass or hypersurface model")
  end

  # Set important fields and return the model
  model.explicit_model_sections = model_sections
  model.defining_classes = defining_classes_of_model
  return model
end



#######################################################
# 4. Constructing models over arbitrary bases
#######################################################

# Construct literature model over arbitrary base
function _construct_literature_model_over_arbitrary_base(model_dict::Dict{String,Any})
  # Construct auxiliary base ring
  @req haskey(model_dict["model_data"], "tunable_sections") "No base coordinates specified for model"
  vars = string.(model_dict["model_data"]["tunable_sections"])
  auxiliary_base_ring, _ = polynomial_ring(QQ, vars, cached=false)

  # Construct the grading of the base ring
  @req haskey(model_dict["model_data"], "classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes") "Database does not specify classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes, but is vital for model construction, so cannot proceed"
  auxiliary_base_grading = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in model_dict["model_data"]["classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes"]]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i:i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  model_dict["model_data"]["classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes"] = auxiliary_base_grading

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
    auxiliary_base_vars = [string(g) for g in symbols(auxiliary_base_ring)]

    # Extract fiber ambient space
    rays = [[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_rays"]]
    max_cones = IncidenceMatrix([[a for a in b] for b in model_dict["model_data"]["fiber_ambient_space_max_cones"]])
    fas = normal_toric_variety(max_cones, rays; non_redundant = true)
    fiber_amb_coordinates = string.(model_dict["model_data"]["fiber_ambient_space_coordinates"])
    set_coordinate_names(fas, fiber_amb_coordinates)

    # Extract the base divisor classes of the fiber coordinates
    fiber_twist_matrix = transpose(matrix(ZZ, (hcat(model_dict["model_data"]["fiber_twist_matrix"]...))))
    @req ncols(fiber_twist_matrix) == length(fiber_amb_coordinates) "Number of fiber coordinate names does not match number of provided fiber gradings"
    @req ncols(fiber_twist_matrix) == n_rays(fas) "Number of rays does not match number of provided fiber gradings"

    # Extract the hypersurface equation
    ambient_ring, _ = polynomial_ring(QQ, vcat(auxiliary_base_vars, fiber_amb_coordinates), cached = false)
    p = eval_poly(model_dict["model_data"]["hypersurface_equation"], ambient_ring)
    
    # Create the model
    model = hypersurface_model(auxiliary_base_vars, auxiliary_base_grading, base_dim, fas, fiber_twist_matrix, p)

  else

    @req false "Model is not a Tate, Weierstrass or hypersurface model"

  end
  return model
end



#######################################################
# 5. Functions for settings attributes of models
#######################################################

function _set_all_attributes(model::AbstractFTheoryModel, model_dict::Dict{String, Any}, model_parameters::Dict{String,<:Any})

  # Metadata
  
  set_attribute!(model, :literature_identifier => model_dict["literature_identifier"])
  
  set_attribute!(model, :model_description => model_dict["model_descriptors"]["description"])
  
  set_attribute!(model, :paper_authors => string.(model_dict["paper_metadata"]["authors"]))
  set_attribute!(model, :paper_buzzwords => string.(model_dict["paper_metadata"]["buzzwords"]))
  set_attribute!(model, :paper_description => model_dict["paper_metadata"]["description"])
  set_attribute!(model, :paper_title => model_dict["paper_metadata"]["title"])

  set_attribute!(model, :arxiv_doi => model_dict["arxiv_data"]["doi"])
  set_attribute!(model, :arxiv_link => model_dict["arxiv_data"]["link"])
  set_attribute!(model, :arxiv_id => model_dict["arxiv_data"]["id"])
  set_attribute!(model, :arxiv_version => model_dict["arxiv_data"]["version"])
  set_attribute!(model, :arxiv_model_equation_number => model_dict["arxiv_data"]["model_location"]["equation"])
  set_attribute!(model, :arxiv_model_page => model_dict["arxiv_data"]["model_location"]["page"])
  set_attribute!(model, :arxiv_model_section => model_dict["arxiv_data"]["model_location"]["section"])
  
  set_attribute!(model, :journal_doi => model_dict["journal_data"]["doi"])
  set_attribute!(model, :journal_link => model_dict["journal_data"]["link"])
  set_attribute!(model, :journal_year => model_dict["journal_data"]["year"])
  set_attribute!(model, :journal_volume => model_dict["journal_data"]["volume"])
  set_attribute!(model, :journal_name => model_dict["journal_data"]["journal"])
  if haskey(model_dict["journal_data"], "report_numbers")
    set_attribute!(model, :journal_report_numbers => string.(model_dict["journal_data"]["report_numbers"]))
  end
  set_attribute!(model, :journal_pages => model_dict["journal_data"]["pages"])
  set_attribute!(model, :journal_model_equation_number => model_dict["journal_data"]["model_location"]["equation"])
  set_attribute!(model, :journal_model_page => model_dict["journal_data"]["model_location"]["page"])
  set_attribute!(model, :journal_model_section => model_dict["journal_data"]["model_location"]["section"])
  
  if haskey(model_dict, "birational_models")
    set_attribute!(model, :birational_literature_models => [str[6:end - 5] for str in model_dict["birational_models"]])
  end
  
  if haskey(model_dict, "associated_models")
    set_attribute!(model, :associated_literature_models => [str[6:end - 5] for str in model_dict["associated_models"]])
  end
  

  # Geometric data

  set_attribute!(model, :partially_resolved, false)

  if haskey(model_dict, "birational_models")
    for m in model_dict["birational_models"]
      model_directory = joinpath(@__DIR__, "Models/")
      model_data = JSON.parsefile(model_directory * m)
      if model_data["model_descriptors"]["type"] == "weierstrass"
        set_attribute!(model, :weierstrass_model => m)
      end
    end
  end

  if haskey(model_dict, "model_parameters")
    set_attribute!(model, :model_parameters => model_parameters)
  end
    
  if haskey(model_dict["model_data"], "resolutions")
    set_attribute!(model, :resolutions => [([string.(c) for c in r[1]], string.(r[2])) for r in model_dict["model_data"]["resolutions"]])
  end

  if haskey(model_dict["model_data"], "weighted_resolutions")
    set_attribute!(model, :weighted_resolutions => [([(string.(c[1]), Int.(c[2])) for c in r[1]], string.(r[2])) for r in model_dict["model_data"]["weighted_resolutions"]])
  end
    
  if haskey(model_dict["model_descriptors"], "global_gauge_group_quotients")
    set_attribute!(model, :global_gauge_group_quotient => map(k -> string.(k), model_dict["model_descriptors"]["global_gauge_group_quotients"]))
  end

  R, _ = polynomial_ring(QQ, collect(keys(explicit_model_sections(model))), cached = false)
  f = hom(R, coordinate_ring(base_space(model)), collect(values(explicit_model_sections(model))))

  if haskey(model_dict["model_data"], "resolution_generating_sections")    
    vs = [[[string.(k) for k in sec] for sec in res] for res in model_dict["model_data"]["resolution_generating_sections"]]
    result = [[[[f(eval_poly(a, R)) for a in b] for b in c] for c in d] for d in vs]
    set_attribute!(model, :resolution_generating_sections => result)
  end
  
  if haskey(model_dict["model_data"], "resolution_zero_sections")
    vs = [[string.(a) for a in b] for b in model_dict["model_data"]["resolution_zero_sections"]]
    result = [[[f(eval_poly(a, R)) for a in b] for b in c] for c in vs]
    set_attribute!(model, :resolution_zero_sections => result)
  end
  
  if haskey(model_dict["model_data"], "weighted_resolution_generating_sections")
    vs = [[[string.(k) for k in sec] for sec in res] for res in model_dict["model_data"]["weighted_resolution_generating_sections"]]
    result = [[[[f(eval_poly(a, R)) for a in b] for b in c] for c in d] for d in vs]
    set_attribute!(model, :weighted_resolution_generating_sections => result)
  end
  
  if haskey(model_dict["model_data"], "weighted_resolution_zero_sections")
    vs = [[string.(a) for a in b] for b in model_dict["model_data"]["weighted_resolution_zero_sections"]]
    result = [[[f(eval_poly(a, R)) for a in b] for b in c] for c in vs]
    set_attribute!(model, :weighted_resolution_zero_sections => result)
  end
  
  if haskey(model_dict["model_data"], "zero_section")
    vs = string.(model_dict["model_data"]["zero_section"])
    set_attribute!(model, :zero_section => [f(eval_poly(l, R)) for l in vs])
  end
  
  if haskey(model_dict["model_data"], "generating_sections")
    vs = map(k -> string.(k), model_dict["model_data"]["generating_sections"])
    set_attribute!(model, :generating_sections => [[f(eval_poly(l, R)) for l in k] for k in vs])
  end
  
  if haskey(model_dict["model_data"], "torsion_sections")
    vs = map(k -> string.(k), model_dict["model_data"]["torsion_sections"])
    set_attribute!(model, :torsion_sections => [[f(eval_poly(l, R)) for l in k] for k in vs])
  end

  if base_space(model) isa NormalToricVariety && ambient_space(model) isa NormalToricVariety
    
    divs = torusinvariant_prime_divisors(ambient_space(model))
    cohomology_ring(ambient_space(model), completeness_check = false)
    cox_gens = symbols(coordinate_ring(ambient_space(model)))

    if haskey(model_dict["model_data"], "zero_section_class") && base_space(model) isa NormalToricVariety
      desired_value = Symbol(string.(model_dict["model_data"]["zero_section_class"]))
      @req desired_value in cox_gens "Specified zero section is invalid"
      index = findfirst(==(desired_value), cox_gens)
      set_attribute!(model, :zero_section_index => index::Int)
      set_attribute!(model, :zero_section_class => cohomology_class(divs[index]))
    end

    if haskey(model_dict["model_data"], "exceptional_classes") && base_space(model) isa NormalToricVariety
      desired_value = string.(model_dict["model_data"]["exceptional_classes"])
      @req issubset(Symbol.(desired_value), cox_gens) "Specified exceptional classes are invalid"
      exceptional_divisor_indices = Vector{Int}()
      for class in desired_value
        index = findfirst(==(Symbol(class)), cox_gens)
        push!(exceptional_divisor_indices, index)
      end
      set_attribute!(model, :exceptional_divisor_indices => exceptional_divisor_indices::Vector{Int})
      set_attribute!(model, :exceptional_classes => [cohomology_class(divs[index]) for index in exceptional_divisor_indices])
    end

  end

  if haskey(model_dict["model_data"], "classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes")
    D = Dict{String, Vector{Int}}()
    tun_sections = model_dict["model_data"]["tunable_sections"]
    M = model_dict["model_data"]["classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes"]
    if typeof(M) != Matrix{Int}
      vars = string.(model_dict["model_data"]["tunable_sections"])
      auxiliary_base_ring, _ = polynomial_ring(QQ, vars, cached=false)
      M = matrix(ZZ, transpose(hcat([[eval_poly(weight, ZZ) for weight in vec] for vec in M]...)))
      M = vcat([[Int(k) for k in M[i:i,:]] for i in 1:nrows(M)]...)
      model_dict["model_data"]["classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes"] = M
    end
    for i in 1:length(tun_sections)
      D[tun_sections[i]] = M[:,i]
    end
    set_attribute!(model, :classes_of_tunable_sections_in_basis_of_Kbar_and_defining_classes, D)
  end
  
  if haskey(model_dict["model_descriptors"], "gauge_algebra")
    algebras = string.(model_dict["model_descriptors"]["gauge_algebra"])
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
    set_attribute!(model, :gauge_algebra => direct_sum(C, LieAlgebra{elem_type(C)}[_construct(g) for g in algebras]))
  end
  
end



#######################################################
# 6. Function to display all known literature models
#######################################################

@doc raw"""
    display_all_literature_models(model_fields::Dict{String,<:Any} = Dict{String,Any}())

Displays all literature models that satisfy the model_fields criteria. The fields currently supported are those occurring in index.json.

```jldoctest; filter = Main.Oscar.doctestfilter_hash_changes_in_1_13()
julia> display_all_literature_models(Dict("gauge_algebra" => ["u(1)", "su(2)", "su(3)"]))
Model 33:
Dict{String, Any}("journal_section" => "3", "arxiv_page" => "67", "arxiv_id" => "1408.4808", "gauge_algebra" => Any["su(3)", "su(2)", "u(1)"], "arxiv_version" => "2", "journal_equation" => "3.141", "journal_page" => "67", "arxiv_equation" => "3.142", "journal_doi" => "10.1007/JHEP01(2015)142", "arxiv_section" => "3", "journal" => "JHEP", "file" => "model1408_4808-11-WSF.json", "arxiv_doi" => "10.48550/arXiv.1408.4808", "model_index" => "33", "type" => "weierstrass")

Model 34:
Dict{String, Any}("journal_section" => "3", "arxiv_page" => "67", "arxiv_id" => "1408.4808", "gauge_algebra" => Any["su(3)", "su(2)", "u(1)"], "arxiv_version" => "2", "journal_equation" => "3.141", "journal_page" => "67", "arxiv_equation" => "3.142", "journal_doi" => "10.1007/JHEP01(2015)142", "arxiv_section" => "3", "journal" => "JHEP", "file" => "model1408_4808-11.json", "arxiv_doi" => "10.48550/arXiv.1408.4808", "model_index" => "34", "type" => "hypersurface")

Model 39:
Dict{String, Any}("journal_section" => "3", "arxiv_page" => "75", "arxiv_id" => "1408.4808", "gauge_algebra" => Any["su(3)", "su(2)", "su(2)", "u(1)"], "arxiv_version" => "2", "journal_equation" => "3.167", "journal_page" => "75", "arxiv_equation" => "3.168", "journal_doi" => "10.1007/JHEP01(2015)142", "arxiv_section" => "3", "journal" => "JHEP", "file" => "model1408_4808-14-WSF.json", "arxiv_doi" => "10.48550/arXiv.1408.4808", "model_index" => "39", "type" => "weierstrass")

Model 40:
Dict{String, Any}("journal_section" => "3", "arxiv_page" => "75", "arxiv_id" => "1408.4808", "gauge_algebra" => Any["su(3)", "su(2)", "su(2)", "u(1)"], "arxiv_version" => "2", "journal_equation" => "3.167", "journal_page" => "75", "arxiv_equation" => "3.168", "journal_doi" => "10.1007/JHEP01(2015)142", "arxiv_section" => "3", "journal" => "JHEP", "file" => "model1408_4808-14.json", "arxiv_doi" => "10.48550/arXiv.1408.4808", "model_index" => "40", "type" => "hypersurface")

Model 45:
Dict{String, Any}("journal_section" => "", "arxiv_page" => "2", "arxiv_id" => "1903.00009", "gauge_algebra" => Any["su(3)", "su(2)", "u(1)"], "arxiv_version" => "3", "journal_equation" => "2", "journal_page" => "2", "arxiv_equation" => "2", "journal_doi" => "10.1103/PhysRevLett.123.101601", "arxiv_section" => "II", "journal" => "Physical Review Letters", "file" => "model1903_00009.json", "arxiv_doi" => "10.48550/arXiv.1903.00009", "model_index" => "45", "type" => "hypersurface")

julia> display_all_literature_models(Dict("gauge_algebra" => "e"))
Model 8:
Dict{String, Any}("journal_section" => "", "arxiv_page" => "49", "arxiv_id" => "1212.2949", "gauge_algebra" => Any["e(6)"], "arxiv_version" => "2", "journal_equation" => "", "journal_page" => "", "arxiv_equation" => "5.1", "journal_doi" => "10.1007/JHEP04(2013)061", "arxiv_section" => "5.1", "journal" => "JHEP", "file" => "model1212_2949-5.json", "arxiv_doi" => "10.48550/arXiv.1212.2949", "model_index" => "8", "type" => "tate")

Model 9:
Dict{String, Any}("journal_section" => "", "arxiv_page" => "49", "arxiv_id" => "1212.2949", "gauge_algebra" => Any["e(7)"], "arxiv_version" => "2", "journal_equation" => "", "journal_page" => "", "arxiv_equation" => "5.7", "journal_doi" => "10.1007/JHEP04(2013)061", "arxiv_section" => "5.1", "journal" => "JHEP", "file" => "model1212_2949-6.json", "arxiv_doi" => "10.48550/arXiv.1212.2949", "model_index" => "9", "type" => "tate")

Model 10:
Dict{String, Any}("journal_section" => "", "arxiv_page" => "49", "arxiv_id" => "1212.2949", "gauge_algebra" => Any["e(8)"], "arxiv_version" => "2", "journal_equation" => "", "journal_page" => "", "arxiv_equation" => "5.13", "journal_doi" => "10.1007/JHEP04(2013)061", "arxiv_section" => "5.1", "journal" => "JHEP", "file" => "model1212_2949-7.json", "arxiv_doi" => "10.48550/arXiv.1212.2949", "model_index" => "10", "type" => "tate")

Model 46:
Dict{String, Any}("journal_section" => "2", "arxiv_page" => "3", "arxiv_id" => "1511.03209", "gauge_algebra" => Any["e(8)", "e(8)", "e(8)", "e(8)", "e(8)", "e(8)", "e(8)", "e(8)", "e(8)", "f(4)", "f(4)", "f(4)", "f(4)", "f(4)", "f(4)", "f(4)", "f(4)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "g(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)", "su(2)"], "arxiv_version" => "3", "journal_equation" => "2.11", "journal_page" => "3", "arxiv_equation" => "2.11", "journal_doi" => "https://doi.org/10.1007/JHEP12(2015)164", "arxiv_section" => "2", "journal" => "JHEP", "file" => "model1511_03209.json", "arxiv_doi" => "10.48550/arXiv.1511.03209", "model_index" => "46", "type" => "tate")
```
"""
function display_all_literature_models(model_fields::Dict{String,<:Any} = Dict{String,Any}())
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  @req issubset(keys(model_fields), keys(file_index[1])) "The inputted criteria aren't supported"
  for field in keys(model_fields)
    if field == "gauge_algebra"
      for s in model_fields["gauge_algebra"]
        filter!(x -> occursin(s, join(x["gauge_algebra"])), file_index)
      end
    else
      filter!(x -> x[field] == model_fields[field], file_index)
    end
  end
  if length(file_index) == 0
    println("No such models found in database")
  end
  sorted_dicts = sort(file_index, by = x -> parse(Int, x["model_index"]))
  for dict in sorted_dicts
    print("Model $(dict["model_index"]):\n")
    print(dict)
    print("\n\n")
  end
end
