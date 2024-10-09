using Aqua

@testset "Aqua.jl" begin
  Aqua.test_all(
    Oscar;
    ambiguities=false,      # TODO: fix ambiguities
    unbound_args=false,     # TODO: fix unbound type parameters
    piracies=false,         # TODO: check the reported methods to be moved upstream
    # Aqua persistent task does not work properly with developed dependencies
    # thus we disable these tests when running in OscarCI:
    persistent_tasks=!haskey(ENV, "oscar_run_tests"),
  )
  @test length(Aqua.detect_unbound_args_recursively(Oscar)) <= 16
end
