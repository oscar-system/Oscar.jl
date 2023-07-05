@doc raw"""
    literature_model(; doi::String="", arxiv_id::String="", version="", equation::String="")

Many models have been created in the F-theory literature.
A significant number of them have even been given specific
names, for instance the "U(1)-restricted SU(5)-GUT model".
This method has access to a database, from which it can
look up such literature models. Currently, you can provide
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
model could be found, or multiple models exist with information matching the
provided information, then an error is raised.

```jldoctest
julia> t = literature_model(arxiv_id = "1109.3454", equation = "3.1")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> v = ambient_space(t)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 0, 0, -2, -3], [0, 0, 0, 1, 0, -2, -3], [0, 0, 0, 0, 1, -2, -3], [0, 1, 0, 0, 0, -2, -3], [0, 0, 1, 0, 0, -2, -3], [0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, -1, -3//2]]

julia> a1,a21,a32,a43,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = blow_up(underlying_toric_variety(v),I)
Normal toric variety

julia> cox_ring(v2)
Multivariate polynomial ring in 9 variables over QQ graded by
  a1 -> [0 0]
  a21 -> [0 0]
  a32 -> [0 0]
  a43 -> [0 0]
  w -> [1 0]
  x -> [1 2]
  y -> [1 3]
  z -> [0 1]
  e -> [-1 0]
```
"""
function literature_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="")
  # Try to find the file with the desired model
  @req any(s -> s != "", [doi, arxiv_id, version, equation]) "No information provided; cannot perform look-up"

  # Read index file
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))

  # Create list of possible candidate files
  candidate_files = Vector{String}()
  for k in 1:length(file_index)
    if all([doi == "" || get(file_index[k], "journal_doi", nothing) == doi,
        arxiv_id == "" || get(file_index[k], "arxiv_id", nothing) == arxiv_id,
        version == "" || get(file_index[k], "arxiv_version", nothing) == version,
        equation == "" || get(file_index[k], "equation", nothing) == equation])
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
      equations = map(d -> get(d["model_location"], "equation", nothing), dicts)
      strings = ["doi: $(dois[i]), arxiv_id: $(ids[i]), version: $(versions[i]), equation: $(equations[i])" for i in 1:length(dicts)]
      "We could not uniquely identify the model. The matched models have the following data:\n$(reduce((s1, s2) -> s1 * "\n" * s2, strings))"
    end)

  # We were able to find a unique model, so read that model
  model_dict = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))

  # Appropriately construct each possible type of model
  if model_dict["model_data"]["type"] == "tate"
    model = _construct_literature_tate_model(model_dict)
  elseif model_dict["model_data"]["type"] == "weierstrass"
    model = _construct_literature_weierstrass_model(model_dict)
  else
    @req false "Model is not a Tate or Weierstrass model"
  end

  # Read in the resolutions and process them to the correct format
  if haskey(model_dict["model_data"], "resolutions")
    resolutions = model_dict["model_data"]["resolutions"]
    resolutions = [[[string.(center) for center in r[1]], string.(r[2])] for r in resolutions]
    set_attribute!(model, :resolutions => resolutions)
  end

  # Remember saved attributes of this model, in particular the resolution sequences known
  if haskey(model_dict["journal_data"], "doi")
    set_attribute!(model, :doi => model_dict["journal_data"]["doi"])
  end
  if haskey(model_dict["arxiv_data"], "id")
    set_attribute!(model, :arxiv_id => model_dict["arxiv_data"]["id"])
  end
  if haskey(model_dict["arxiv_data"], "version")
    set_attribute!(model, :version => model_dict["arxiv_data"]["version"])
  end
  if haskey(model_dict["model_location"], "equation")
    set_attribute!(model, :equation_number => model_dict["model_location"]["equation"])
  end
  if haskey(model_dict["model_data"], "description")
    set_attribute!(model, :description => model_dict["model_data"]["description"])
  end
  if haskey(model_dict["arxiv_data"], "link")
    set_attribute!(model, :link => model_dict["arxiv_data"]["link"])
  end

  # Return the model
  return model
end

# Constructs Tate model from given Tate literature model
function _construct_literature_tate_model(model_dict::Dict{String,Any})
  auxiliary_base_ring, _ = PolynomialRing(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)

  base_dim = get(model_dict["model_data"], "base_dim", 3)

  a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_base_ring)
  a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_base_ring)
  a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_base_ring)
  a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_base_ring)
  a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_base_ring)

  return global_tate_model([a1, a2, a3, a4, a6], auxiliary_base_ring, base_dim)
end

# Constructs Weierstrass model from given Weierstrass literature model
function _construct_literature_weierstrass_model(model_dict::Dict{String,Any})
  auxiliary_base_ring, _ = PolynomialRing(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)

  base_dim = get(model_dict["model_data"], "base_dim", 3)

  f = eval_poly(get(model_dict["model_data"], "f", "0"), auxiliary_base_ring)
  g = eval_poly(get(model_dict["model_data"], "g", "0"), auxiliary_base_ring)

  return weierstrass_model(f, g, auxiliary_base_ring, base_dim)
end