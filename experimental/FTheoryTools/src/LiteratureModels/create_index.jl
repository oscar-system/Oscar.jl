# Function to automatically generate index file from existing .json files
# named model*.json in the Models directory

# Currently, this stores the following data for quick searching of models:
#    * model_index
#    * arXiv data
#        -arXiv ID
#        -arXiv DOI
#        -arXiv version number
#    * publication data
#        -arXiv DOI
#        -journal name
#    * paper meta_data
#        -author list
#        -paper title
#    * model location within paper
#        -section number that introduces model
#        -equation number that introduces model
#        -page number that introduces model (in arXiv version of paper, currently)
#    * model_descriptors
#        -type

# internal helper
function _model_indices()
  return JSON.parsefile(joinpath(@__DIR__, "model_indices.json"))
end

function _create_literature_model_index()
  model_directory = joinpath(@__DIR__, "Models")
  models = readdir(model_directory)
  filter!(startswith("model"), models)

  index = Vector{Dict{String,Union{String,Vector{Any}}}}()
  model_indices = _model_indices()
  for model in models
    model_data = JSON.parsefile(joinpath(model_directory, model))
    model_index = get(model_indices, model, "")
    if model_index == ""
      model_index = string.(maximum([parse(Int, x) for x in values(model_indices)]) + 1)
      model_indices[model] = model_index
    end

    model_index_dict = Dict("model_index" => model_index)

    arxiv_data = get(model_data, "arxiv_data", false)
    if arxiv_data != false
      arxiv_dict = Dict(
        "arxiv_id" => get(arxiv_data, "id", ""),
        "arxiv_doi" => get(arxiv_data, "doi", ""),
        "arxiv_version" => get(arxiv_data, "version", ""),
        "arxiv_section" => get(arxiv_data["model_location"], "section", ""),
        "arxiv_equation" => get(arxiv_data["model_location"], "equation", ""),
        "arxiv_page" => get(arxiv_data["model_location"], "page", ""))
    else
      arxiv_dict = Dict{String,String}()
    end

    journal_data = get(model_data, "journal_data", false)
    if journal_data != false
      journal_dict = Dict(
        "journal_doi" => get(journal_data, "doi", ""),
        "journal" => get(journal_data, "journal", ""),
        "journal_section" => get(journal_data["model_location"], "section", ""),
        "journal_equation" => get(journal_data["model_location"], "equation", ""),
        "journal_page" => get(journal_data["model_location"], "page", ""))
    else
      journal_dict = Dict{String,String}()
    end

    meta_data = get(model_data, "paper_meta_data", false)
    if meta_data != false
      meta_data_dict = Dict(
        "authors" => get(meta_data, "authors", [""]), "title" => get(meta_data, "title", "")
      )
    else
      meta_data_dict = Dict{String,Union{String,Vector{Any}}}()
    end

    model_descriptor_data = get(model_data, "model_descriptors", false)
    if model_descriptor_data != false
      model_descriptor_data_dict = Dict(
        "type" => get(model_descriptor_data, "type", [""]),
        "gauge_algebra" => get(model_descriptor_data, "gauge_algebra", [""]),
      )
    else
      model_descriptor_data_dict = Dict{String,Union{String,Vector{Any}}}()
    end

    index_entry = merge(
      model_index_dict, arxiv_dict, journal_dict, meta_data_dict, model_descriptor_data_dict
    )
    index_entry["file"] = model
    if !isempty(index_entry)
      push!(index, index_entry)
    end
  end

  #Check if any models have been removed and update the index accordingly
  if issetequal(collect(keys(model_indices)), models) == false
    difference = setdiff(collect(keys(model_indices)), models)
    for key in difference
      delete!(model_indices, key)
    end
  end

  open(joinpath(@__DIR__, "index.json"), "w") do file
    JSON.print(file, index)
  end
  open(joinpath(@__DIR__, "model_indices.json"), "w") do file
    JSON.print(file, model_indices, 2)
  end
end
