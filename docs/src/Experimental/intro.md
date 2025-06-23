# Adding new projects to experimental

## Purpose
The folder `experimental` is for code that is candidate for being added to
Oscar. In particular, this includes the following cases:
- Code from an external package that should be moved to Oscar
- Implementing a new feature from scratch
- In general: code whose API has not stabilized yet

The code in `experimental` is supposed to be mathematically correct,
`experimental` is a staging area for new features whose interface is still to
be stabilized. Also code is allowed to reside in `experimental` while it is
brought up to Oscar standard.

!!! danger "Dependencies"
    - Code from `src` must never use code from `experimental`
    - Say there are two projects `A` and `B` in `experimental`, and `B` depends
      on `A`. That means that `B` cannot be moved to `src` before `A`. Worse:
      If `A` gets abandoned, `B` might share that fate. So please consider
      carefully in such situations.

## Structure
For an example of the structure for a new project in `experimental` have a look
at project folders, i.e. `experimental/PROJECT_NAME`, that have subfolders
`docs`, `src`, and `test`. The general structure is
```
experimental/PROJECT_NAME/
├── README.md
├── docs
│   ├── doc.main
│   └── src
│       └── DOCUMENTATION.md
├── src
│   └── PROJECT_NAME.jl
└── test
    └── *.jl
```
The file `src/PROJECT_NAME.jl` and at least one `.jl` file in the `test/`
directory are mandatory and are used by Oscar.jl to find your code and tests.
If there is a `test/runtests.jl` then only this file is executed during
testing, otherwise all `.jl` files will be run automatically (in a random
order).

The file `docs/doc.main` is used for integrating your documentation in the
Oscar manual under the `Experimental` section. Optionally please provide a
`README.md` describing your project and its goals, especially if you are
starting from scratch and don't have any documentation yet.

### Useful functions for development
Apart from the hints in the [Introduction for new developers](@ref), there are some more specialized functions for the structure of the `experimental` folder.
```@docs
Oscar.test_experimental_module
```

## Procedure for adding a new feature
Ideally we envision the procedure to follow along the following lines.

1. The new feature is implemented in the `experimental` folder.
2. For external authors, a maintainer is assigned to guide the authors such
   that the implementation adheres to the [Developer Style Guide](@ref) and the
   [Design Decisions](@ref).
   Please get in touch with us as soon as possible, preferably on the [OSCAR
   Slack](https://oscar.computeralgebra.de/community/#slack).
3. The new feature is tested thoroughly, other people are invited to test the
   new feature.
4. In the end there are three possibilities:
   1. The feature is considered done and moved into `src` as is.
   2. The feature is discarded, e.g., because it cannot be maintained.
   3. Parts of the feature are moved into `src`, others are discarded.

## Criteria for acceptance

The main criteria for acceptance are:
1. The code adheres to the [Developer Style Guide](@ref) and the [Design
   Decisions](@ref).
2. The new code is well tested.
3. It is clear who maintains the new code, i.e. the original authors commit to
   maintaining the code in the future.

