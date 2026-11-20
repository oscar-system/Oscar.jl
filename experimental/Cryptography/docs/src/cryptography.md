# Post-quantum cryptography

This section describes implementations of cryptographic schemes that are
believed to be secure against attacks by quantum computers: the multivariate
signature scheme UOV, the lattice-based encryption scheme LWE, several
hash-based signature schemes and code-based encryption schemes. They are meant
for experiments and teaching. The implementations are neither constant-time nor
hardened against side channels and must not be used to protect real data.

## The common interface

All schemes use the same functions, which dispatch on the parameter object `p`
of the scheme, the first argument:

- `keygen(p)` returns a key pair `(pk, sk)`,
- `sign(p, msg, sk)` signs a message and `verify(p, sig, msg, pk)` checks a
  signature (signature schemes), `open(p, sm, pk)` returns the message of a signed
  message `sm` (a message followed by its signature) or `nothing`,
- `encrypt(p, pk, msg)` and `decrypt(p, sk, ct)` (encryption schemes).

Parameter sets are created with `instantiate_uov`, `instantiate_lwe`,
`instantiate_hbs` and `instantiate_mceliece`, or taken from the predefined ones
such as `uov_1p`, `lwe_toy`, `hbs_mss_toy` or `mceliece_toy`.

## UOV signatures

The Unbalanced Oil and Vinegar signature scheme with the parameter sets of the
NIST submission, see [BCDGKPPS.EA25](@cite).

```@docs
instantiate_uov
keygen(::UOV)
sign(::UOV, ::Vector{UInt8}, ::Vector{UInt8})
verify(::UOV, ::Vector{UInt8}, ::Vector{UInt8}, ::Vector{UInt8})
open(::UOV, ::Vector{UInt8}, ::Vector{UInt8})
```

## Learning with errors

Public key encryption based on the learning with errors problem, see
[Reg05](@cite).

```@docs
instantiate_lwe
keygen(::LWE)
encrypt(::LWE, ::Any, ::AbstractVector{UInt8})
decrypt(::LWE, ::Any, ::Any)
```

## Hash-based signatures

Signature schemes whose security rests on the properties of a hash function only
(here SHAKE256): the one-time signatures of Lamport [Lam79](@cite) and Winternitz
(WOTS+, [Hls13](@cite)), Merkle trees [Mer90](@cite), the few-time signatures
HORST and FORS, and their combination into signature schemes with a single Merkle
tree or with a hypertree in the style of SPHINCS [BHHLNPSS.EA15](@cite) and
SPHINCS+ [BHKNRS19](@cite). The parameter sets are not compatible with the NIST
standard FIPS 205.

### Signatures

```@docs
instantiate_hbs
keygen(::HBS)
sign(::HBS, ::Vector{UInt8}, ::Vector{UInt8})
verify(::HBS, ::Vector{UInt8}, ::Vector{UInt8}, ::Vector{UInt8})
open(::HBS, ::Vector{UInt8}, ::Vector{UInt8})
```

### One-time signatures

```@docs
lamport_keygen
lamport_sign
lamport_verify
```

### Merkle trees

```@docs
merkle_tree
merkle_root
merkle_auth_path
merkle_verify
```

## Code-based encryption

Public key encryption whose security rests on the hardness of decoding random
linear codes. The secret key is a binary Goppa code, whose errors can be
corrected with Patterson's algorithm [Pat75](@cite). The parameter sets are not
compatible with the NIST submission Classic McEliece.

### McEliece

The scheme of McEliece [McE78](@cite): a message block is encoded with a
scrambled generator matrix of the code and random errors are added.

```@docs
instantiate_mceliece
keygen(::McEliece)
encrypt(::McEliece, ::FqMatrix, ::AbstractVector{UInt8})
decrypt(::McEliece, ::Any, ::FqMatrix)
```

## Hash function

```@docs
shake_256
```
