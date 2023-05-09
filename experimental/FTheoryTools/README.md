# Welcome to FTheoryTools

## Goal

We aim to automate a number of recurring and (at least in part) tedious computations in F-theory model building
as described in large detail in [Wei18](@cite). Specifically, we focus on the following setups:
* 4d F-theory compactifications,
* defined by a global and singular Weierstrass model as codimension 1 locus of a toric ambient space $Y$,
* which can be crepantly resolved.
Some of the techniques/algorithms extend naturally to more general settings. For example, it is not at all necessary to restrict to 4-dimensional settings (or alternatively, base spaces of dimension 3). Indeed, the current implementation does allow for *arbitrary* base dimension. For more extensions that we might address in the future, please take a look at the section "possible future extensions" below.

We aim for the following workflow:
* User input:
    * Weierstrass polynomial $P_W$,
    * Data defining the toric ambient space $Y$ (if applicable),
    * Choice of resolved phase (if applicable),
    * Generating sections (for $\operatorname{U}(1)$ symmetries).
* Output:
    * Singular loci in codimension 1, 2, and 3,
    * Defining data of resolved geometry,
    * (Pictures of) fibre diagrams of resolved fibre over the originally singular loci, including intersections of $\operatorname{U}(1)$-sections,
    * Gauge group,
    * Topological data (e.g., Euler number).


## Status

This project just began, and is therefore in its experimental stage. Upcoming tasks include, but are not limited, to the following:
* The irrelevant ideal, SR ideal, and ideal of linear relations may need to be modified when the exceptional coordinate "e" is included in the blowup. They are currently set up to work when e is eliminated.
* Decide whether to stick with global blowups or use charts.
* Consolidate notation about sections and line bundles in the documentation.
* The Kodaira type function assumes that the singular locus is given by a single coordinate.
* The Kodaira type function only works for codimension 1.
* Modify the _ambient_space_from_base function to return ring maps from the base to the ambient space and all other constructors/types to appropriately carry that around, then fix the corresponding bit of code in analyze_fibers.


## Tutorial

We provide a [tutorial for FTheoryTools in OSCAR](https://nbviewer.org/github/oscar-system/oscar-website/blob/gh-pages/tutorials/FTheoryToolsInOSCAR.ipynb).


## Possible future extensions

Future extensions include, but are not necessarily limited to, the following:
* Specify a $G_4$-flux and work out the chiral spectra,
* Specify a gauge potential and work out (candidates for) the line bundles whose cohomologies encode the vector-like spectra,
* Other singularity types (non-minimal, terminal, etc.,
* Base blowups for singularity resolution.


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Martin Bies](https://martinbies.github.io/),
* [Andrew Turner](https://apturner.net/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on GitHub](https://github.com/oscar-system/Oscar.jl).


## Acknowledgements

We appreciate insightful discussions with [Mirjam Cvetič](https://live-sas-physics.pantheon.sas.upenn.edu/people/standing-faculty/mirjam-cvetic). The work of Andrew Turner is supported by DOE (HEP) Award DE-SC001352.
