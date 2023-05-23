# global seed value for oscar to allow creating deterministic random number
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

