@doc raw"""
    literature_tate_model(; arxiv_id::String="", equ_nr::String="", description::String="")

Many Tate models have been created in the F-theory literature.
A significant number of them have even been given specific
names, for instance the "U(1)-restricted SU(5)-GUT model".
This method has access to a database, from which it can
look up such literature models. Currently, you can provide
any combination of the following three optional arguments
to the method `literature_tate_model`:
* `doi`: A string representing the DOI of the publication that
introduced the model in question.
* `equ_nr`: A string representing the number of the equation, which introduced
the model in question.
* `description`: A string with a description of the model in question.
For papers, that were published on the arxiv, we can instead of the `doi` also
provide the following:
* `arxiv_id`: For papers published on the arxiv, one can also provide a string
which represents the arxiv-identifier of the paper which introduced the model
in question.
* `version`: A string representing the version of the arxiv-upload.
The method `literature_tate_model` attempts to find a model in our data base
for which the provided data matches the information in our record. If no such
model could be found, or multiple models exist with information matching the
provided information, then the following error is raised: "We could not uniquely
 identify the model".

```jldoctest
julia> t = literature_tate_model(arxiv_id = "1109.3454", equ_nr = "3.5")
Global Tate model over a not fully specified base -- SU(5)xU(1) restricted Tate model based on arxiv paper 1109.3454 (equ. 3.5)

julia> v = ambient_space(t)
Scheme of a toric variety with fan spanned by RayVector{QQFieldElem}[[1, 0, 0, 0, 0, 0, -2, -3], [0, 0, 0, 0, 1, 0, -2, -3], [0, 0, 0, 0, 0, 1, -2, -3], [0, 1, 0, 0, 0, 0, -2, -3], [0, 0, 1, 0, 0, 0, -2, -3], [0, 0, 0, 1, 0, 0, -2, -3], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, -1, -3//2]]

julia> a10,a21,a32,a43,a65,w,x,y,z = gens(cox_ring(v));

julia> I = ideal([x,y,w]);

julia> v2 = blow_up(underlying_toric_variety(v),I)
Normal toric variety

julia> cox_ring(v2)
Multivariate Polynomial Ring in 10 variables a10, a21, a32, a43, ..., e over Rational Field graded by
  a10 -> [0 0]
  a21 -> [0 0]
  a32 -> [0 0]
  a43 -> [0 0]
  a65 -> [0 0]
  w -> [1 0]
  x -> [1 2]
  y -> [1 3]
  z -> [0 1]
  e -> [-1 0]
```
"""
function literature_tate_model(; doi::String="", arxiv_id::String="", version::String="", equ_nr::String="", description::String="")
  
  # try to find the file with the desired model
  @req (arxiv_id != "" || equ_nr != "" || description != "") "No information provided -- cannot perform look-up"
    
  # read index file
  file_index = JSON.parsefile(joinpath(@__DIR__, "index.json"))
  
  # create list of possible candidate files
  candidate_files = Vector{String}()
  for k in 1:length(file_index)
    if doi != "" && haskey(file_index[k], "doi") && file_index[k]["doi"] != doi
      continue
    end
    if arxiv_id != "" && haskey(file_index[k], "arxiv_id") && file_index[k]["arxiv_id"] != arxiv_id
      continue
    end
    if version != "" && haskey(file_index[k], "version") && file_index[k]["version"] != version
      continue
    end
    if equ_nr != "" && haskey(file_index[k], "equ_nr") && file_index[k]["equ_nr"] != equ_nr
      continue
    end
    if description != "" && haskey(file_index[k], "description") && file_index[k]["description"] != description
      continue
    end
    push!(candidate_files, string(file_index[k]["file"]))
  end
  
  # check if we found exactly one file, i.e. we were able to identify the model uniquely
  @req length(candidate_files) == 1 "We could not uniquely identify the model"
  
  # we were able to find a unique model, so read that model
  model = JSON.parsefile(joinpath(@__DIR__, "Models/" * candidate_files[1]))
  
  # Construct the model in question
  auxiliary_base_ring, vars = PolynomialRing(QQ, string.(model["variables"]), cached=false)
  a1 = zero(auxiliary_base_ring)
  if haskey(model, "a1") && length(model["a1"]) > 0
    a1 = prod([vars[p[1]]^p[2] for p in model["a1"]])
  end
  a2 = zero(auxiliary_base_ring)
  if haskey(model, "a2") && length(model["a2"]) > 0
    a2 = prod([vars[p[1]]^p[2] for p in model["a2"]])
  end
  a3 = zero(auxiliary_base_ring)
  if haskey(model, "a3") && length(model["a3"]) > 0
    a3 = prod([vars[p[1]]^p[2] for p in model["a3"]])
  end
  a4 = zero(auxiliary_base_ring)
  if haskey(model, "a4") && length(model["a4"]) > 0
    a4 = prod([vars[p[1]]^p[2] for p in model["a4"]])
  end
  a6 = zero(auxiliary_base_ring)
  if haskey(model, "a6") && length(model["a6"]) > 0
    a6 = prod([vars[p[1]]^p[2] for p in model["a6"]])
  end
  t = global_tate_model([a1, a2, a3, a4, a6], auxiliary_base_ring, Int.(model["dim"]))
  
  # Read in the resolutions and process them to the correct format
  if haskey(model, "resolutions")
    resolutions = model["resolutions"];
    resolutions = [[string.(x) for x in r] for r in resolutions]
    set_attribute!(t, :resolutions => resolutions)
  end
  
  # Remember saved attributes of this model, in particular the resolution sequences known
  if haskey(model, "doi")
    set_attribute!(t, :doi => model["doi"])
  end
  if haskey(model, "arxiv_id")
    set_attribute!(t, :arxiv_id => model["arxiv_id"])
  end
  if haskey(model, "version")
    set_attribute!(t, :version => model["version"])
  end
  if haskey(model, "equ_nr")
    set_attribute!(t, :equ_nr => model["equ_nr"])
  end
  if haskey(model, "description")
    set_attribute!(t, :description => model["description"])
  end
  if haskey(model, "link")
    set_attribute!(t, :link => model["link"])
  end
  
  # Return the model
  return t
end
