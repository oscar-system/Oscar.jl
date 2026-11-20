# ---------------------------------------------------------------------------
# Code-based encryption schemes, whose security rests on the hardness of
# decoding random linear codes: the schemes of McEliece and Niederreiter with
# binary Goppa codes.
#
# The implementation is meant for experiments and teaching. It is neither
# constant-time nor protected against chosen-ciphertext attacks.
# ---------------------------------------------------------------------------
using Random: AbstractRNG, RandomDevice, shuffle

include("goppa.jl")
