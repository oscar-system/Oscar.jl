function create_literature_model_index()
  model_directory = joinpath(@__DIR__, "Models/")
  models = readdir(model_directory)
  filter!(s -> startswith(s, "model"), models)

  index = Vector{Dict{String,Union{String,Vector{Any}}}}()
  for model in models
    model_data = JSON.parsefile(model_directory * model)

    arxiv_data = get(model_data, "arxiv_data", false)
    if arxiv_data !== false
      arxiv_dict = Dict("arxiv_id" => get(arxiv_data, "id", ""), "arxiv_doi" => get(arxiv_data, "doi", ""), "arxiv_version" => get(arxiv_data, "version", ""))
    else
      arxiv_dict = Dict{String,String}()
    end

    journal_data = get(model_data, "journal_data", false)
    if journal_data !== false
      journal_dict = Dict("journal_doi" => get(journal_data, "doi", ""), "journal" => get(journal_data, "journal", ""))
    else
      journal_dict = Dict{String,String}()
    end

    metadata = get(model_data, "paper_metadata", false)
    if metadata !== false
      metadata_dict = Dict("authors" => get(metadata, "authors", [""]), "title" => get(metadata, "title", ""))
    else
      metadata_dict = Dict{String,Union{String,Vector{Any}}}()
    end

    location_data = get(model_data, "model_location", false)
    if location_data !== false
      location_dict = Dict("section" => get(location_data, "section", ""), "equation" => get(location_data, "equation", ""), "page" => get(location_data, "page", ""))
    else
      location_dict = Dict{String,String}()
    end

    index_entry = merge(arxiv_dict, journal_dict, metadata_dict, location_dict)
    if !isempty(index_entry)
      push!(index, index_entry)
    end
  end

  open("index.json", "w") do file
    JSON.print(file, index)
  end
end
