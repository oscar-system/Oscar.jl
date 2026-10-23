using Aqua

@testset "Aqua.jl" begin
  # stale_deps and persistent_tasks load Oscar in a new process, which is
  # slow with code coverage enabled; the other CI jobs run them
  coverage = Base.JLOptions().code_coverage != 0
  Aqua.test_all(
    Oscar;
    ambiguities=false,      # TODO: fix ambiguities
    piracies=false,         # TODO: check the reported methods to be moved upstream
    stale_deps = !coverage && (; ignore = [:Compat]),  # conditionally loaded
    # Aqua persistent task does not work properly with developed dependencies
    # thus we disable these tests when running in OscarCI:
    persistent_tasks=!coverage && !haskey(ENV, "oscar_run_tests"),
  )
end
