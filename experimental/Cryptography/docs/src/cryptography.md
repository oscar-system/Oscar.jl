# Post-quantum cryptography

This section describes implementations of cryptographic schemes that are
believed to be secure against attacks by quantum computers: the multivariate
signature scheme UOV, the lattice-based encryption scheme LWE and several
hash-based signature schemes. They are meant for experiments and teaching. The
implementations are neither constant-time nor hardened against side channels and
must not be used to protect real data.

## UOV signatures

The Unbalanced Oil and Vinegar signature scheme with the parameter sets of the
NIST submission, see [BCDGKPPS.EA25](@cite).

```@docs
instantiate_uov
uov_keygen
uov_sign
uov_verify
uov_open
```

## Learning with errors

Public key encryption based on the learning with errors problem, see
[Reg05](@cite).

```@docs
instantiate_lwe
lwe_keygen
lwe_encrypt
lwe_decrypt
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
hbs_keygen
hbs_sign
hbs_verify
hbs_open
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

## Hash function

```@docs
shake_256
```
