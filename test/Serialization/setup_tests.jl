using JSONSchema, Oscar.JSON
# This code is needed for multiple test files that may end up on different workers.
# Thus, this needs to be conditionally included in each of these test files.

# we only need to define this once

if !isdefined(Main, :mrdi_schema)
  schemajson = JSON.parsefile(joinpath(Oscar.oscardir, "data", "schema.json"))
  # replace remote ref to polymake schema with local copy to avoid network access
  if schemajson["\$defs"]["data"]["oneOf"][5]["\$ref"] == "https://polymake.org/schemas/data.json"
    # this needs to be an absolute path
    schemajson["\$defs"]["data"]["oneOf"][5]["\$ref"] = "file://$(joinpath(Oscar.oscardir,"data","polymake.json"))"
  else
    error("mardi schema: please update json-path for polymake schema reference")
  end
  mrdi_schema = Schema(schemajson)
end

if !isdefined(Main, :test_save_load_roundtrip) || isinteractive()
  function test_save_load_roundtrip(func, path, original::T;
                                    params=nothing, check_func=nothing,
                                    serializer::Oscar.Serialization.OscarSerializer=Oscar.Serialization.JSONSerializer(),
                                    kw...) where {T}
    is_multi_file = serializer isa Oscar.Serialization.MultiFileSerializer
    is_multi_file_ref_serializer = serializer isa Oscar.Serialization.MultiFileRefSerializer

    # save and load from a file
    filename = is_multi_file_ref_serializer ? joinpath(path, "original") : joinpath(path, "original.json")
    save(filename, original; serializer=serializer, kw...)
    loaded = load(filename; params=params, serializer=serializer, kw...)

    @test loaded isa T
    func(loaded)

    if !is_multi_file
      # save and load from a file without saving references
      no_refs = Oscar.Serialization.JSONSerializer(serialize_refs=false)
      save(filename, original; serializer=no_refs, kw...)
      @test !any(line -> contains(line, "_refs"), eachline(filename))
      loaded = load(filename; params=params, serializer=no_refs, kw...)

      @test loaded isa T
      func(loaded)

      # save and load from an IO buffer
      io = IOBuffer()
      save(io, original; serializer=serializer, kw...)
      seekstart(io)
      loaded = load(io; params=params, serializer=serializer, kw...)

      @test loaded isa T
      func(loaded)

      # save and load from an IO buffer, with prescribed type
      io = IOBuffer()
      save(io, original; serializer=serializer, kw...)
      seekstart(io)
      loaded = load(io; type=T, params=params, serializer=serializer, kw...)

      @test loaded isa T
      func(loaded)
    end

    # test loading on a empty state
    save(filename, original; serializer=serializer, kw...)
    Oscar.reset_global_serializer_state()
    loaded = load(filename; params=params, serializer=serializer, kw...)
    @test loaded isa T

    # test passing TypeAndParams
    save(filename, original; kw...)
    Oscar.reset_global_serializer_state()
    loaded = load(filename; params=Oscar.type_and_params(original), kw...)
    @test loaded isa T

    # test schema
    schema_file = is_multi_file_ref_serializer ? filename * ".mrdi" : filename
    jsondict = JSON.parsefile(schema_file)
    @test validate(mrdi_schema, jsondict) == nothing

    isnothing(check_func) || @test check_func(loaded)
  end
end
