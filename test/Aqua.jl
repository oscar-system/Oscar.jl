using Aqua

@testset "Aqua.jl" begin
  Aqua.test_all(
    Oscar;
    ambiguities=false,      # TODO: fix ambiguities
    piracies=false,         # TODO: check the reported methods to be moved upstream
    stale_deps = (; ignore = [:Compat]),  # conditionally loaded
    # Aqua persistent task does not work properly with developed dependencies
    # thus we disable these tests when running in OscarCI:
    persistent_tasks=!haskey(ENV, "oscar_run_tests"),
  )
end
