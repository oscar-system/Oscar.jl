if !isdefined(Main, :test_1_4_0_upgrade) || isinteractive()
  function test_1_4_0_upgrade(;exclude::Vector{String}=String[], only::Union{Vector{String},Nothing}=nothing)
    artifact_toml = Oscar.LazyArtifacts.find_artifacts_toml(Oscar.oscardir)
    Oscar.LazyArtifacts.ensure_artifact_installed("version-1-3-0-files", artifact_toml)
    _hash = Oscar.LazyArtifacts.artifact_hash("version-1-3-0-files", artifact_toml)
    dir = Oscar.LazyArtifacts.artifact_path(_hash)
    type_folders = joinpath(dir, "version_1_3_0_files")
    types_done = Set{String}()
    for dir_name in readdir(type_folders)
      type_str = split(dir_name, "-")[1]
      type_str in exclude && continue
      !isnothing(only) && !(type_str in only) && continue
      push!(types_done, type_str)
      T = Oscar.reverse_type_map[type_str]
      @testset "$T" begin
        @testset "$filename" for filename in readdir(joinpath(type_folders, dir_name))
          @test load(joinpath(type_folders, dir_name, filename)) isa T
        end
      end
    end
    if !isnothing(only) && !issubset(only, types_done)
      error("The following types are not found in version_1_3_0_files: $(setdiff(only, types_done))")
    end
  end
end
