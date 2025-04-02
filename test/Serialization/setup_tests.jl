using JSONSchema, Oscar.JSON
# This code is needed for multiple test files that may end up on different workers.
# Thus, this needs to be conditionally included in each of these test files.

# we only need to define this once

if !isdefined(Main, :mrdi_schema)
  schemajson = JSON.parsefile(joinpath(Oscar.oscardir, "data", "schema.json"))
  # replace remote ref to polymake schema with local copy to avoid network access
  if schemajson["\$defs"]["data"]["oneOf"][4]["\$ref"] == "https://polymake.org/schemas/data.json"
    # this needs to be an absolute path
    schemajson["\$defs"]["data"]["oneOf"][4]["\$ref"] = "file://$(joinpath(Oscar.oscardir,"data","polymake.json"))"
  else
    error("mardi schema: please update json-path for polymake schema reference")
  end
  mrdi_schema = Schema(schemajson)
end

if !isdefined(Main, :test_save_load_roundtrip) || isinteractive()
  function test_save_load_roundtrip(func, path, original::T;
                                    params=nothing, check_func=nothing, kw...) where {T}
    # save and load from a file
    filename = joinpath(path, "original.json")
    save(filename, original; kw...)
    loaded = load(filename; params=params, kw...)
    
    @test loaded isa T
    func(loaded)

    # save and load from an IO buffer
    io = IOBuffer()
    save(io, original; kw...)
    seekstart(io)
    loaded = load(io; params=params, kw...)

    @test loaded isa T
    func(loaded)

    # save and load from an IO buffer, with prescribed type
    io = IOBuffer()
    save(io, original; kw...)
    seekstart(io)
    loaded = load(io; type=T, params=params, kw...)

    @test loaded isa T
    func(loaded)

    # test loading on a empty state
    save(filename, original; kw...)
    Oscar.reset_global_serializer_state()
    loaded = load(filename; params=params, kw...)
    @test loaded isa T

    # test schema
    jsondict = JSON.parsefile(filename)
    @test validate(mrdi_schema, jsondict) == nothing

    isnothing(check_func) || @test check_func(loaded)
  end
end
