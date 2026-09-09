```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Extra-long tests

## Why do we have extra-long tests?

The test suite that runs for every push and pull request should provide broad
coverage while completing in a reasonable amount of time. Tests involving
particularly large examples or expensive computations should therefore not be
added to the regular continuous integration test suite.

At the same time, some computations are important enough that their results
should be checked regularly even though they take too long for the regular
continuous integration test suite. Examples include large regression tests,
computations supporting mathematical publications, and tests ensuring that
published results remain reproducible with current versions of OSCAR. For such
computations, we provide the extra-long test suite, which is run once per day.


## Adding an extra-long test

To add tests to the extra-long test suite, place the relevant `@testset` (or
`@testset`s) in a dedicated Julia test file, add this file to the `test` directory
just like any other OSCAR test file, and register its path in the
`:extra_long` entry of `test_subsets` in `test/runtests.jl`.


## Running the tests on a pull request

To run the extra-long test suite for a pull request, add the `extra-long` label
to the pull request and trigger the workflow.  This can be done, for example, by
pushing a new commit or by closing and reopening the pull request.
