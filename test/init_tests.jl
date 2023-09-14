# This file is dedicated for code that needs to be run on every worker before
# running the tests. This is useful for test-code that is called from multiple files
# that may run on different workers.

function test_save_load_roundtrip(func, path, original::T; params=nothing) where T
  # save and load from a file
  filename = joinpath(path, "original.json")
  save(filename, original)
  loaded = load(filename; params=params)
  if T <: Vector
    @test loaded isa Vector
  else
    @test loaded isa T
  end
  func(loaded)

  # save and load from an IO buffer
  io = IOBuffer()
  save(io, original)
  seekstart(io)
  loaded = load(io; params=params)

  if T <: Vector
    @test loaded isa Vector
  else
    @test loaded isa T
  end
  func(loaded)

  # save and load from an IO buffer, with prescribed type
  io = IOBuffer()
  save(io, original)
  seekstart(io)
  loaded = load(io, type=T, params=params)
  if T <: Vector
    @test loaded isa Vector
  else
    @test loaded isa T
  end
  func(loaded)

  # test loading on a empty state
  save(filename, original)
  reset_global_serializer_state()
  loaded = load(filename; params=params)
end
