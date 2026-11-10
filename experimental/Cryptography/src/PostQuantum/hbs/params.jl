# ---------------------------------------------------------------------------
# Parameter sets. A parameter set describes either a scheme with a single
# Merkle tree over Winternitz keys (`:mss`, stateful) or a scheme with a
# hypertree of several layers of Merkle trees, whose leaves sign the roots of
# the layer below, and whose bottom layer signs the public keys of few-time
# signatures (`:hypertree`, stateless). The few-time signature is HORST or FORS.
# ---------------------------------------------------------------------------

struct HBS
    scheme::Symbol
    fts::Symbol
    n::Int
    w::Int
    h::Int
    d::Int
    hp::Int
    len::Int
    tau::Int
    k::Int
    x::Int
    pk_sz::Int
    sk_sz::Int
    sig_sz::Int
    name::String
end

function Base.show(io::IO, p::HBS)
    print(io, "Hash-based signature parameters ", p.name, " (", p.scheme)
    p.fts == :none || print(io, ", ", p.fts)
    print(io, ", n = ", p.n, ", w = ", p.w, ", h = ", p.h)
    p.d == 1 || print(io, ", d = ", p.d)
    print(io, ")")
end

# Number of bytes of a signature of the few-time signature scheme.
function _fts_sig_size(fts::Symbol, n::Int, tau::Int, k::Int, x::Int)
    fts == :fors && return k * (n + tau * n)
    return k * (n + (tau - x) * n) + 2^x * n
end

@doc raw"""
    instantiate_hbs(n::Int, w::Int, h::Int, name::String; d::Int = 1, fts::Symbol = :none, tau::Int = 0, k::Int = 0, x::Int = 0) -> HBS

Create a parameter set for a hash-based signature scheme with hashes of `n`
bytes and Winternitz parameter `w`. There are two kinds of schemes:
- With `fts = :none` (the default) the scheme uses a single Merkle tree of
  height `h`, whose leaves are the public keys of $2^h$ Winternitz keys. It can
  sign $2^h$ messages and is stateful: signing changes the secret key.
- With `fts = :horst` or `fts = :fors` the scheme uses a hypertree of total
  height `h` with `d` layers of Merkle trees, so that $2^h$ few-time keys can
  sign. A message is signed by a few-time key chosen by its digest, which is
  then certified by the layers of the hypertree. For HORST, `tau` is the height
  of the tree, `k` the number of revealed values and `x` the number of top
  levels of the tree that are part of the signature. For FORS, there are `k`
  trees of height `tau`, and `x` is not used. The scheme is stateless.

The parameters must satisfy the following restrictions:
- `n >= 1` and `w` is a power of two between 2 and 256;
- `h >= 1` and `d` divides `h`; for `fts = :none`, we need `d = 1` and `h <= 24`;
- for `fts != :none`, `tau`, `k >= 1`, `tau <= 24`, and for HORST `0 <= x <= tau`.

The label `name` is only used for printing. The predefined sets `hbs_mss_toy`,
`hbs_mss_light`, `hbs_horst_toy`, `hbs_horst_light`, `hbs_fors_toy`,
`hbs_fors_light` and `hbs_fors_standard` are created this way.

# Examples
```jldoctest
julia> instantiate_hbs(16, 16, 10, "my-mss")
Hash-based signature parameters my-mss (mss, n = 16, w = 16, h = 10)

julia> instantiate_hbs(16, 16, 12, "my-fors"; d = 3, fts = :fors, tau = 8, k = 16)
Hash-based signature parameters my-fors (hypertree, fors, n = 16, w = 16, h = 12, d = 3)

julia> instantiate_hbs(16, 15, 10, "bad")
ERROR: ArgumentError: w must be a power of two between 2 and 256, got 15
[...]
```
"""
function instantiate_hbs(n::Int, w::Int, h::Int, name::String; d::Int = 1, fts::Symbol = :none, tau::Int = 0, k::Int = 0, x::Int = 0)
    @req n >= 1 "n must be positive, got $n"
    @req 2 <= w <= 256 && ispow2(w) "w must be a power of two between 2 and 256, got $w"
    @req h >= 1 "h must be positive, got $h"
    @req fts in (:none, :horst, :fors) "fts must be :none, :horst or :fors, got :$fts"
    len = _wots_len(n, w)
    if fts == :none
        @req d == 1 "a scheme with one Merkle tree needs d = 1, got $d"
        @req h <= 24 "the height of a single tree must be at most 24, got $h"
        @req tau == 0 && k == 0 && x == 0 "tau, k and x are only used with a few-time signature"
        sig_sz = 8 + n + len * n + h * n
        return HBS(:mss, fts, n, w, h, 1, h, len, 0, 0, 0, 2n, 4n + 8, sig_sz, name)
    end
    @req d >= 1 && h % d == 0 "d must divide h, got h = $h and d = $d"
    @req 1 <= tau <= 24 && k >= 1 "tau must be between 1 and 24 and k must be positive, got tau = $tau and k = $k"
    @req h <= 62 "the total height must be at most 62, got $h"
    if fts == :horst
        @req 0 <= x <= tau "x must be between 0 and tau, got $x"
    else
        @req x == 0 "x is not used with FORS"
    end
    hp = h ÷ d
    @req hp <= 24 "the height h/d of one tree must be at most 24, got $hp"
    sig_sz = n + _fts_sig_size(fts, n, tau, k, x) + d * (len * n + hp * n)
    return HBS(:hypertree, fts, n, w, h, d, hp, len, tau, k, x, 2n, 4n, sig_sz, name)
end

# Parameter sets. They are illustrative and are not claimed to provide any
# concrete security level. The toy sets are tiny and fast, `hbs_fors_standard`
# is modelled on the set SPHINCS+-128f of [BHKNRS19](@cite) with total height 60
# instead of 66 (indices must fit into an `Int`) and is not compatible with it.
const hbs_mss_toy = instantiate_hbs(8, 4, 4, "hbs-mss-toy")
const hbs_mss_light = instantiate_hbs(16, 16, 10, "hbs-mss-light")
const hbs_horst_toy = instantiate_hbs(8, 4, 6, "hbs-horst-toy"; d = 2, fts = :horst, tau = 6, k = 8, x = 2)
const hbs_horst_light = instantiate_hbs(16, 16, 12, "hbs-horst-light"; d = 3, fts = :horst, tau = 10, k = 24, x = 4)
const hbs_fors_toy = instantiate_hbs(8, 4, 6, "hbs-fors-toy"; d = 2, fts = :fors, tau = 4, k = 8)
const hbs_fors_light = instantiate_hbs(16, 16, 12, "hbs-fors-light"; d = 3, fts = :fors, tau = 8, k = 16)
const hbs_fors_standard = instantiate_hbs(16, 16, 60, "hbs-fors-standard"; d = 20, fts = :fors, tau = 6, k = 33)

const hbs_all = [hbs_mss_toy, hbs_mss_light, hbs_horst_toy, hbs_horst_light, hbs_fors_toy, hbs_fors_light, hbs_fors_standard]
