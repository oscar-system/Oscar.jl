# global seed value for OSCAR to allow creating deterministic random number
# generators
# we initialize this here with a random value as well to allow use during
# precompilation
const rng_seed = Ref{UInt32}(rand(UInt32))

@doc raw"""
    get_seed()

Return the current random seed that is used for calls to `Oscar.get_seeded_rng`.
"""
get_seed() = return rng_seed[]

@doc raw"""
    set_seed!(s::Integer)

Set a new global seed for all subsequent calls to `Oscar.get_seeded_rng`.
"""
function set_seed!(s::Integer)
  rng_seed[] = convert(UInt32, s)
end

@doc raw"""
    get_seeded_rng()

Return a new random number generator object of type MersenneTwister which is
seeded with the global seed value `Oscar.rng_seed`. This can be used for
the testsuite to get more stable output and running times. Using a separate RNG
object for each testset (or file) makes sure previous uses don't affect later
testcases. It could also be useful for some randomized algorithms.
The seed will be initialized with a random value during OSCAR startup but can
be set to a fixed value with `Oscar.set_seed!` as it is done in `runtests.jl`.
"""
get_seeded_rng() = return MersenneTwister([get_seed()])


"""
    randseed!(s::Integer)

Reset the random seed for many global RNGs used by OSCAR. Currently that includes:
- the default Julia random source, `Random.default_rng()`
- the default FLINT random source
- the default GAP random source `GlobalMersenneTwister`
- the legacy GAP random source `GlobalRandomSource`
- the Singular random sources

Note that no special steps are taken for thread-local random sources (e.g. the
default Julia random source currently is even task local). Code using multiple
threads or multiple tasks may need to call this function from each thread or
task.

This function is *not* supposed to be called from regular code in Oscar.jl.
Rather it is meant to help with debugging issues that depend on random
sources. E.g. if a computation only fails sometimes, depending on randomness,
it may be possible to make it reproducible by a trick like this:
```julia
  seed=0
  while true
    Oscar.randseed!(seed)
    run_tricky_computation()
    seed += 1
  end
```
Once this throws an exception, we can simply record the value of `seed` at
that time, and use that to make the error reproducible:
```julia
  Oscar.randseed!(1234)
  run_tricky_computation()  # now always throws an exception
```

# Examples

```jldoctest
julia> Oscar.randseed!(1234)

julia> rand(1:1000)
219

julia> GAP.Globals.Random(1,1000)
163
```
"""
randseed!(s::Integer) = randseed!(UInt32(s))

function randseed!(s::UInt32)
  # reset Julia random source on the current thread
  Random.seed!(s)

  # reset FLINT random source on the current thread
  Nemo.randseed!(s)

  # reset GAP random sources
  GAP.randseed!(s)

  # reset Singular random sources
  Singular.randseed!(s)

  return nothing
end
