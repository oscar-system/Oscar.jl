```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Performance tools

The OSCAR timing dashboard records and visualizes timing data of
existing OSCAR tests as benchmarks. It can be used to identify performance regressions,
inspect benchmark runs, and compare timing results across OSCAR commits.

The dashboard is available at <https://speed.oscar-system.org>.

The source code for the dashboard is maintained at [oscar-system/oscar-timing](https://github.com/oscar-system/oscar-timing).

A dedicated server periodically checks the OSCAR `master` branch for new commits and
benchmarks them in chronological order. Thus, new timing data is generated only
when commits have been added to `master`. If no commits are added for several days,
the dashboard data remains unchanged and may be marked as stale; this does not by
itself indicate a problem with the benchmarking infrastructure.
