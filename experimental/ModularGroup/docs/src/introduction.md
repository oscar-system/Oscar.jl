```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

This package provides methods for computing with finite-index subgroups of the modular groups
${\rm SL}_2(\mathbb{Z})$ and ${\rm PSL}_2(\mathbb{Z})$.

A convenient way to represent finite-index subgroups of ${\rm SL}_2(\mathbb{Z})$ is by specifying the action
of generator matrices of ${\rm SL}_2(\mathbb{Z})$ on the right cosets by right multiplication. For example,
one could choose the generators

```math
S = \begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix} \qquad T = \begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}
```

and represent a subgroup as a tuple of transitive permutations $(\sigma_S, \sigma_T)$ describing the action
of $S$ and $T$. This is exactly the way this package internally treats such subgroups. We use the convention
that $1$ corresponds to the coset of the identity matrix. Note that such a representation as a tuple of
permutations is only unique up to relabeling of the cosets, i.e. up to simultaneous conjugation (fixing the
$1$ coset by our convention).

## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what
this means. See also the dedicated `README.md` for details.

## Contact

Please direct questions about this part of OSCAR to the following people:

- [Sebastian Engelhardt](mailto:seen00001@stud.uni-saarland.de)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on GitHub](https://www.oscar-system.org/community/#how-to-report-issues).
