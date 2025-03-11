if !isdefined(Main, :test_1_4_0_upgrade) || isinteractive()
  function test_1_4_0_upgrade(;
    exclude::Vector#={<:Union{String,Pair{String,AbstractVector{Int}}}}=#=String[],
    only::Union{Vector#={<:Union{String,Pair{String,AbstractVector{Int}}}}=#,Nothing}=nothing,
  )  
    exclude_types = [type for type in exclude if type isa String]
    exclude_type_ids = Dict{String, Vector{Int}}(first(pair) => collect(last(pair)) for pair in exclude if pair isa Pair)
    
    if !isnothing(only)
      only_types = [type for type in only if type isa String]
      only_type_ids = Dict{String, Vector{Int}}(first(pair) => collect(last(pair)) for pair in only if pair isa Pair)
    end

    artifact_toml = Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)
    Oscar.LazyArtifacts.ensure_artifact_installed("version-1-3-0-files", artifact_toml)
    _hash = Oscar.LazyArtifacts.artifact_hash("version-1-3-0-files", artifact_toml)
    dir = Oscar.LazyArtifacts.artifact_path(_hash)
    type_folders = joinpath(dir, "version_1_3_0_files")
    for dir_name in readdir(type_folders)
      type_str = split(dir_name, "-")[1]
      type_str in exclude_types && continue
      !isnothing(only) && !(type_str in only_types) && !(haskey(only_type_ids, type_str)) && continue
      T = Oscar.reverse_type_map[type_str]
      @testset "$T" begin
        for filename in readdir(joinpath(type_folders, dir_name))
          if haskey(exclude_type_ids, type_str)
            id = parse(Int, split(filename, ".")[1])
            id in exclude_type_ids[type_str] && continue
          end
          if !isnothing(only) && haskey(only_type_ids, type_str)
            id = parse(Int, split(filename, ".")[1])
            !(id in only_type_ids[type_str]) && continue
          end
          @testset "$filename" begin
            @test load(joinpath(type_folders, dir_name, filename)) isa T
          end
        end
      end
    end
    Oscar.reset_global_serializer_state()
  end
end
