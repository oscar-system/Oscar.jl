# Implementing features from scratch

For new features whose interface naturally will change a lot during the
implementation phase, Oscar.jl has a staging area, the `experimental` folder.

1. The new feature is implemented in the `experimental` folder.
2. Additionally to the original authors, a maintainer is assigned to guide the
   authors such that the implementation adheres to the [Developer Style
   Guide](@ref) and the [Design Decisions](@ref).
3. The new feature is tested thoroughly, other people are invited to test the
   new feature.
4. In the end there are three possibilities:
   1. The feature is considered done and moved into `src` as is.
   2. The feature is discarded, since it cannot be maintained.
   3. Parts of the feature are moved into `src`, others are discarded.

## Criteria for acceptance

The main criteria for acceptance are:
1. The code adheres to the [Developer Style Guide](@ref) and the [Design
   Decisions](@ref).
2. The new code is well tested.
3. It is clear who maintains the new code, i.e. the original authors commit to
   maintaining the code in the future.

