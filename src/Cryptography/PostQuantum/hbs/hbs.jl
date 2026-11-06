# ---------------------------------------------------------------------------
# Hash-based signature schemes: one-time signatures (Lamport, Winternitz),
# Merkle trees, few-time signatures (HORST, FORS) and the combinations of
# these into signature schemes with a single Merkle tree or with a hypertree
# of several Merkle trees. All hashing is done with SHAKE256.
#
# The implementation is meant for experiments and teaching. It is neither
# constant-time nor compatible with the NIST standard FIPS 205.
# ---------------------------------------------------------------------------
using Random: AbstractRNG, RandomDevice

include("hash.jl")
include("merkle.jl")
include("lamport.jl")
include("wots.jl")
