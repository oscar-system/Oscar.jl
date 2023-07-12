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
Assuming that the first row of the given grading is the grading under Kbar

Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arXiv paper 1109.3454 Eq. (3.1)

julia> v = ambient_space(t)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 2, -2, 0, 1], [-1, -3//2, 1//2, 0, 0], [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, -1, 1//3], [0, 0, 0, 1, -1//2], [0, 0, 0, 0, 1]]

julia> a1,a21,a32,a43,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = blow_up(underlying_toric_variety(v),I)
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

  # We were able to find a unique model, so read that model)
  model_dict = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))

  # Appropriately construct each possible type of model
  if model_dict["model_data"]["type"] == "tate"
    model = _construct_literature_tate_model(model_dict)
  elseif model_dict["model_data"]["type"] == "weierstrass"
    model = _construct_literature_weierstrass_model(model_dict)
  else
    @req false "Model is not a Tate or Weierstrass model"
  end

  # Remember saved attributes of this model
  base_ring = cox_ring(base_space(model)) # THIS CURRENTLY ASSUMES THE BASE IS TORIC, SHOULD FIX
  if haskey(model_dict["arxiv_data"], "id")
    set_attribute!(model, :arxiv_id => model_dict["arxiv_data"]["id"])
  end
  if haskey(model_dict["arxiv_data"], "doi")
    set_attribute!(model, :arxiv_doi => model_dict["arxiv_data"]["doi"])
  end
  if haskey(model_dict["arxiv_data"], "link")
    set_attribute!(model, :arxiv_link => model_dict["arxiv_data"]["link"])
  end
  if haskey(model_dict["arxiv_data"]["model_location"], "equation")
    set_attribute!(model, :arxiv_model_equation_number => model_dict["arxiv_data"]["model_location"]["equation"])
  end
  if haskey(model_dict["arxiv_data"]["model_location"], "page")
    set_attribute!(model, :arxiv_model_page => model_dict["arxiv_data"]["model_location"]["page"])
  end
  if haskey(model_dict["arxiv_data"]["model_location"], "section")
    set_attribute!(model, :arxiv_model_section => model_dict["arxiv_data"]["model_location"]["section"])
  end
  if haskey(model_dict["arxiv_data"], "version")
    set_attribute!(model, :arxiv_version => model_dict["arxiv_data"]["version"])
  end
  if haskey(model_dict, "associated_models")
    # This removes 'model' from the front and '.json' from the end of the file name
    set_attribute!(model, :associated_literature_models => [str[6:end - 5] for str in model_dict["associated_models"]])
  end
  if haskey(model_dict["model_data"], "generating_sections")
    set_attribute!(model, :generating_sections => [[eval_poly(coord, base_ring) for coord in gen_sec] for gen_sec in model_dict["model_data"]["generating_sections"]])
  end
  if haskey(model_dict["journal_data"], "doi")
    set_attribute!(model, :journal_doi => model_dict["journal_data"]["doi"])
  end
  if haskey(model_dict["journal_data"], "link")
    set_attribute!(model, :journal_link => model_dict["journal_data"]["link"])
  end
  if haskey(model_dict["journal_data"]["model_location"], "equation")
    set_attribute!(model, :journal_model_equation_number => model_dict["journal_data"]["model_location"]["equation"])
  end
  if haskey(model_dict["journal_data"]["model_location"], "page")
    set_attribute!(model, :journal_model_page => model_dict["journal_data"]["model_location"]["page"])
  end
  if haskey(model_dict["journal_data"]["model_location"], "section")
    set_attribute!(model, :journal_model_section => model_dict["journal_data"]["model_location"]["section"])
  end
  if haskey(model_dict["journal_data"], "pages")
    set_attribute!(model, :journal_pages => model_dict["journal_data"]["pages"])
  end
  if haskey(model_dict["journal_data"], "report_numbers")
    set_attribute!(model, :journal_report_numbers => string.(model_dict["journal_data"]["report_numbers"]))
  end
  if haskey(model_dict["journal_data"], "volume")
    set_attribute!(model, :journal_volume => model_dict["journal_data"]["volume"])
  end
  if haskey(model_dict["journal_data"], "year")
    set_attribute!(model, :journal_year => model_dict["journal_data"]["year"])
  end
  # This removes 'model' from the front and '.json' from the end of the file name
  set_attribute!(model, :literature_identifier => candidate_files[1][6:end - 5])
  if haskey(model_dict["model_data"], "description")
    set_attribute!(model, :model_description => model_dict["model_data"]["description"])
  end
  if haskey(model_dict["paper_metadata"], "authors")
    set_attribute!(model, :paper_authors => string.(model_dict["paper_metadata"]["authors"]))
  end
  if haskey(model_dict["paper_metadata"], "buzzwords")
    set_attribute!(model, :paper_buzzwords => string.(model_dict["paper_metadata"]["buzzwords"]))
  end
  if haskey(model_dict["paper_metadata"], "description")
    set_attribute!(model, :paper_description => model_dict["paper_metadata"]["description"])
  end
  if haskey(model_dict["paper_metadata"], "title")
    set_attribute!(model, :paper_title => model_dict["paper_metadata"]["title"])
  end
  if haskey(model_dict, "related_models")
    # This removes 'model' from the front and '.json' from the end of the file name
    set_attribute!(model, :related_literature_models => [str[6:end - 5] for str in model_dict["related_models"]])
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
  if haskey(model_dict["model_data"], "zero_section")
    set_attribute!(model, :zero_section => [eval_poly(coord, base_ring) for coord in model_dict["model_data"]["zero_section"]])
  end
  # Return the model
  return model
end

# Constructs Tate model from given Tate literature model
function _construct_literature_tate_model(model_dict::Dict{String,Any})
  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = PolynomialRing(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)
  
  base_dim = get(model_dict["model_data"], "base_dim", 3)
  
  a1 = eval_poly(get(model_dict["model_data"], "a1", "0"), auxiliary_base_ring)
  a2 = eval_poly(get(model_dict["model_data"], "a2", "0"), auxiliary_base_ring)
  a3 = eval_poly(get(model_dict["model_data"], "a3", "0"), auxiliary_base_ring)
  a4 = eval_poly(get(model_dict["model_data"], "a4", "0"), auxiliary_base_ring)
  a6 = eval_poly(get(model_dict["model_data"], "a6", "0"), auxiliary_base_ring)
  
  @req haskey(model_dict["model_data"], "auxiliary_base_grading") "Currently, only literature models over arbitrary bases are supported"
  auxiliary_base_grading = matrix(ZZ, transpose(hcat(model_dict["model_data"]["auxiliary_base_grading"]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  
  return global_tate_model(auxiliary_base_ring, auxiliary_base_grading, base_dim, [a1, a2, a3, a4, a6])
end

# Constructs Weierstrass model from given Weierstrass literature model
function _construct_literature_weierstrass_model(model_dict::Dict{String,Any})
  @req haskey(model_dict["model_data"], "base_coordinates") "No base coordinates specified for model"
  auxiliary_base_ring, _ = PolynomialRing(QQ, string.(model_dict["model_data"]["base_coordinates"]), cached=false)
  
  base_dim = get(model_dict["model_data"], "base_dim", 3)
  
  f = eval_poly(get(model_dict["model_data"], "f", "0"), auxiliary_base_ring)
  g = eval_poly(get(model_dict["model_data"], "g", "0"), auxiliary_base_ring)
  
  @req haskey(model_dict["model_data"], "auxiliary_base_grading") "Currently, only literature models over arbitrary bases are supported"
  auxiliary_base_grading = matrix(ZZ, transpose(hcat(model_dict["model_data"]["auxiliary_base_grading"]...)))
  auxiliary_base_grading = vcat([[Int(k) for k in auxiliary_base_grading[i,:]] for i in 1:nrows(auxiliary_base_grading)]...)
  
  return weierstrass_model(auxiliary_base_ring, auxiliary_base_grading, base_dim, f, g)
end
