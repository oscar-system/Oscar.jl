using JSONSchema, Oscar.JSON
# This code is needed for multiple test files that may end up on different workers.
# Thus, this needs to be conditionally included in each of these test files.

# we only need to define this once

if !isdefined(Main, :test_save_load_roundtrip)
  mrdi_schema = Schema(JSON.parsefile(joinpath(Oscar.oscardir, "data", "schema.json")))

  function test_save_load_roundtrip(func, path, original::T; params=nothing) where {T}
    # save and load from a file
    filename = joinpath(path, "original.json")
    save(filename, original)
    loaded = load(filename; params=params)
    
    @test loaded isa T
    func(loaded)

    # save and load from an IO buffer
    io = IOBuffer()
    save(io, original)
    seekstart(io)
    loaded = load(io; params=params)

    @test loaded isa T
    func(loaded)

    # save and load from an IO buffer, with prescribed type
    io = IOBuffer()
    save(io, original)
    seekstart(io)
    loaded = load(io; type=T, params=params)

    @test loaded isa T
    func(loaded)

    # test loading on a empty state
    save(filename, original)
    Oscar.reset_global_serializer_state()
    loaded = load(filename; params=params)
    @test loaded isa T

    # test schema
    jsondict = JSON.parsefile(filename)
    @test validate(mrdi_schema, jsondict) == nothing
  end

end
