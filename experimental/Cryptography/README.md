# Cryptography

## Aims

This package collects implementations of cryptographic schemes for experiments
and teaching, in particular post-quantum schemes:

- the multivariate signature scheme UOV (`src/PostQuantum/uov`),
- the lattice-based encryption scheme LWE (`src/PostQuantum/lwe`),
- hash-based signatures (`src/PostQuantum/hbs`): Lamport and Winternitz one-time
  signatures, Merkle trees, the few-time signatures HORST and FORS, and
  signature schemes with a single Merkle tree or a hypertree,
- the hash function SHAKE256 (`src/keccak.jl`) that the hash-based schemes use.

## Status

The implementations are neither constant-time nor hardened against side channels
and are not meant to protect real data. The hash-based schemes are not
compatible with the NIST standard FIPS 205.
