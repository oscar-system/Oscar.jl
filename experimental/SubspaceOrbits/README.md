# Orbits on subspaces

## Aims

Compute representatives, lengths and stabilizers of the orbits of a matrix
group `G` over a prime field `F` on the `k`-dimensional subspaces of `F^d`,
for cases with many subspaces.

A subspace is identified with its basis matrix in reduced row echelon form,
which is encoded as an integer between 1 and the number of `k`-subspaces;
see O'Brien, *The p-group generation algorithm*, J. Symbolic Comput. 9
(1990), Sections 3.3-3.5. The set of subspaces already found is then a
bitmap with one bit per subspace. Besides that the orbit being computed is
stored with 4 or 8 bytes per subspace, and for stabilizers a parent edge of
16 bytes per block.

- `src/labels.jl`: the bijection between subspaces and integers.
- `src/label_actions.jl`: the image of an integer under a matrix. Rows of
  basis matrices are machine words, with one bit per entry for `p = 2` and
  several bits per entry for odd `p`.
- `src/orbits.jl`: the orbit algorithm. For a solvable normal subgroup `N`
  of `G` (the solvable radical) the `N`-orbits are computed along a
  polycyclic generating sequence, which takes one image per subspace. The
  `N`-orbits are blocks for `G`, and the remaining generators of `G` are
  applied to whole blocks. The order in which the points of an orbit are
  found encodes a transversal, from which stabilizers are read off.
- `src/interface.jl`: the functions for OSCAR matrix groups.

## Status

The exported functions are `subspace_orbit_representatives_and_lengths` and
`subspace_orbit_representatives_and_stabilizers`. The latter returns the
same kind of result as `orbit_representatives_and_stabilizers` and is meant
to become a further `algorithm` of that function.

Limits:

- Prime fields only, with `p < 2^16`.
- At most `2^36` subspaces of the requested dimension.
- GAP computes with `G` through its action on vectors, so nothing is asked
  of GAP if `F^d` has more than `2^24` vectors: lengths are then computed by
  the plain orbit algorithm with the given generators, and stabilizers are
  refused.
- For lengths alone the order and the solvable radical of `G` are computed
  only if the radical can save at least `2^27` images, since GAP may need
  far longer for them than the plain orbit algorithm for the orbits.
- The code for the image of a subspace is generated for each pair `(k, d)`;
  the first call for a new pair compiles for a few seconds.

Not done:

- Non-prime fields.
- Stabilizers if `F^d` has more than `2^24` vectors. They need the order of
  `G` and membership tests in `G`; `RecognizeGroup` from the GAP package
  recog may become a way to get these without the action on vectors.
- For a group whose orbits consist of many blocks, stabilizers are generated
  by Schreier generators of random edges until the known order is reached.
  Each of them needs a scan of the orbit and a membership test in a
  permutation group, and their number is not bounded in advance.
