@doc raw"""
    literature_tate_model(; doi::String="", arxiv_id::String="", version="", equation::String="")

Many Tate models have been created in the F-theory literature.
A significant number of them have even been given specific
names, for instance the "U(1)-restricted SU(5)-GUT model".
This method has access to a database, from which it can
look up such literature models. Currently, you can provide
any combination of the following three optional arguments
to the method `literature_tate_model`:
* `doi`: A string representing the DOI of the publication that
introduced the model in question.
* `equation`: A string representing the number of the equation that introduced
the model in question.
For papers, that were published on the arxiv, we can instead of the `doi` also
provide the following:
* `arxiv_id`: For papers published on the arxiv, one can also provide a string
that represents the arxiv identifier of the paper that introduced the model
in question.
* `version`: A string representing the version of the arxiv upload.
The method `literature_tate_model` attempts to find a model in our database
for which the provided data matches the information in our record. If no such
model could be found, or multiple models exist with information matching the
provided information, then the following error is raised: "We could not uniquely
 identify the model".

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equation = "3.1")
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
function literature_tate_model(; doi::String="", arxiv_id::String="", version::String="", equation::String="")
  
  # try to find the file with the desired model
  @req (arxiv_id != "" || equation != "" || description != "") "No information provided; cannot perform look-up"
    
  # read index file
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  
  # create list of possible candidate files
  candidate_files = Vector{String}()
  for k in 1:length(file_index)
    if doi != "" && haskey(file_index[k], "journal_doi") && file_index[k]["journal_doi"] != doi
      continue
    end
    if arxiv_id != "" && haskey(file_index[k], "arxiv_id") && file_index[k]["arxiv_id"] != arxiv_id
      continue
    end
    if version != "" && haskey(file_index[k], "arxiv_version") && file_index[k]["arxiv_version"] != version
      continue
    end
    if equation != "" && haskey(file_index[k], "equation") && file_index[k]["equation"] != equation
      continue
    end
    push!(candidate_files, string(file_index[k]["file"]))
  end
  
  # check if we found exactly one file, i.e. we were able to identify the model uniquely
  @req length(candidate_files) == 1 "We could not uniquely identify the model"
  
  # we were able to find a unique model, so read that model
  model = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))

  # ensure that the found model is a Tate model
  @req model["model_data"]["type"] == "tate" "The given model is not a Tate model"
  
  # Construct the model in question
  auxiliary_base_ring, _ = PolynomialRing(QQ, string.(model["model_data"]["base_coordinates"]), cached=false)
  a1 = zero(auxiliary_base_ring)
  if haskey(model["model_data"], "a1")
    a1 = eval_poly(model["model_data"]["a1"], auxiliary_base_ring)
  end
  a2 = zero(auxiliary_base_ring)
  if haskey(model["model_data"], "a2")
    a2 = eval_poly(model["model_data"]["a2"], auxiliary_base_ring)
  end
  a3 = zero(auxiliary_base_ring)
  if haskey(model["model_data"], "a3")
    a3 = eval_poly(model["model_data"]["a3"], auxiliary_base_ring)
  end
  a4 = zero(auxiliary_base_ring)
  if haskey(model["model_data"], "a4")
    a4 = eval_poly(model["model_data"]["a4"], auxiliary_base_ring)
  end
  a6 = zero(auxiliary_base_ring)
  if haskey(model["model_data"], "a6")
    a6 = eval_poly(model["model_data"]["a6"], auxiliary_base_ring)
  end
  t = global_tate_model([a1, a2, a3, a4, a6], auxiliary_base_ring, Int.(model["model_data"]["base_dim"]))
  
  # Read in the resolutions and process them to the correct format
  if haskey(model["model_data"], "resolutions")
    resolutions = model["model_data"]["resolutions"];
    resolutions = [[[string.(center) for center in r[1]], string.(r[2])] for r in resolutions]
    set_attribute!(t, :resolutions => resolutions)
  end
  
  # Remember saved attributes of this model, in particular the resolution sequences known
  if haskey(model["journal_data"], "doi")
    set_attribute!(t, :doi => model["journal_data"]["doi"])
  end
  if haskey(model["arxiv_data"], "id")
    set_attribute!(t, :arxiv_id => model["arxiv_data"]["id"])
  end
  if haskey(model["arxiv_data"], "version")
    set_attribute!(t, :version => model["arxiv_data"]["version"])
  end
  if haskey(model["model_location"], "equation")
    set_attribute!(t, :equation_number => model["model_location"]["equation"])
  end
  if haskey(model["model_data"], "description")
    set_attribute!(t, :description => model["model_data"]["description"])
  end
  if haskey(model["arxiv_data"], "link")
    set_attribute!(t, :link => model["arxiv_data"]["link"])
  end
  
  # Return the model
  return t
end
