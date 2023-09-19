# Function to automatically generate index file from existing .json files
# named model*.json in the Models directory

# Currently, this stores the following data for quick searching of models:
#    * arXiv data
#        -arXiv ID
#        -arXiv DOI
#        -arXiv version number
#    * publication data
#        -arXiv DOI
#        -journal name
#    * paper metadata
#        -author list
#        -paper title
#    * model location within paper
#        -section number that introduces model
#        -equation number that introduces model
#        -page number that introduces model (in arXiv version of paper, currently)


function _create_literature_model_index()
  model_directory = joinpath(@__DIR__, "Models/")
  models = readdir(model_directory)
  filter!(s -> startswith(s, "model"), models)

  index = Vector{Dict{String,Union{String,Vector{Any}}}}()
  for model in models
    model_data = JSON.parsefile(model_directory * model)

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

    metadata = get(model_data, "paper_metadata", false)
    if metadata != false
      metadata_dict = Dict("authors" => get(metadata, "authors", [""]), "title" => get(metadata, "title", ""))
    else
      metadata_dict = Dict{String,Union{String,Vector{Any}}}()
    end

    index_entry = merge(arxiv_dict, journal_dict, metadata_dict)
    index_entry["file"] = model
    if !isempty(index_entry)
      push!(index, index_entry)
    end
  end

  open(joinpath(@__DIR__,"index.json"), "w") do file
    JSON.print(file, index)
  end
end
