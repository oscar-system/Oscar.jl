# ---------------------------------------------------------------------------
# Code-based encryption schemes, whose security rests on the hardness of
# decoding random linear codes: the schemes of McEliece and Niederreiter with
# binary Goppa codes and the scheme HQC with quasi-cyclic codes.
#
# The implementation is meant for experiments and teaching. It is neither
# constant-time nor protected against chosen-ciphertext attacks.
# ---------------------------------------------------------------------------
using Random: AbstractRNG, RandomDevice, randperm, shuffle

include("bits.jl")
include("goppa.jl")
include("mceliece.jl")
include("niederreiter.jl")
include("hqc_codes.jl")
