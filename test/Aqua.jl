using Aqua

@testset "Aqua.jl" begin
  Aqua.test_all(
    Oscar;
    ambiguities=false,      # TODO: fix ambiguities
    unbound_args=false,     # TODO: fix unbound type parameters
    piracies=false,         # TODO: check the reported methods to be moved upstream
  )
  @test length(Aqua.detect_unbound_args_recursively(Oscar)) <= 16
end
