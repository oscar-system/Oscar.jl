import Downloads
import CodecZlib
import Tar

if !isdefined(Main, :serialization_upgrade_test_path) ||
  !isdir(Main.serialization_upgrade_test_path) ||
  !isfile(joinpath(Main.serialization_upgrade_test_path, "LICENSE.md"))

  serialization_upgrade_test_path = let commit_hash = "d3d6a7876c135ce7d3945679597ebe7100e87a58"
    tarball = Downloads.download("https://github.com/oscar-system/serialization-upgrade-tests/archive/$(commit_hash).tar.gz")

    destpath = open(CodecZlib.GzipDecompressorStream, tarball) do io
      Tar.extract(io)
    end
    joinpath(destpath, "serialization-upgrade-tests-$(commit_hash)")
  end
end

if !isdefined(Main, :test_upgrade_folder) || isinteractive()
  function test_upgrade_folder(folder::String;
    exclude::Vector#={<:Union{String,Pair{String,AbstractVector{Int}}}}=#=String[],
    only::Union{Vector#={<:Union{String,Pair{String,AbstractVector{Int}}}}=#,Nothing}=nothing,
  )  
    exclude_types = [type for type in exclude if type isa String]
    exclude_type_ids = Dict{String, Vector{Int}}(first(pair) => collect(last(pair)) for pair in exclude if pair isa Pair)
    
    if !isnothing(only)
      only_types = [type for type in only if type isa String]
      only_type_ids = Dict{String, Vector{Int}}(first(pair) => collect(last(pair)) for pair in only if pair isa Pair)
    end

    type_folders = joinpath(Main.serialization_upgrade_test_path, folder)
    for dir_name in readdir(type_folders)
      type_str = split(dir_name, "-")[1]
      type_str in exclude_types && continue
      !isnothing(only) && !(type_str in only_types) && !(haskey(only_type_ids, type_str)) && continue
      T = Oscar.Serialization.reverse_type_map[type_str]
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
