using Aqua

@testset "Aqua.jl" begin
  # stale_deps and persistent_tasks load Oscar in a new process, which
  # takes minutes in the Julia 1.10 CI jobs; the other CI jobs run them
  slow_reload = VERSION < v"1.11"
  Aqua.test_all(
    Oscar;
    ambiguities=false,      # TODO: fix ambiguities
    piracies=false,         # TODO: check the reported methods to be moved upstream
    stale_deps = !slow_reload && (; ignore = [:Compat]),  # conditionally loaded
    # Aqua persistent task does not work properly with developed dependencies
    # thus we disable these tests when running in OscarCI:
    persistent_tasks=!slow_reload && !haskey(ENV, "oscar_run_tests"),
  )
end
