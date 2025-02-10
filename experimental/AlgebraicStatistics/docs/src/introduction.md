```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

Algebraic Statistics uses tools from algebra, geometry, and combinatorics to solve problems in statistics. In particular, since many parametric statistical models are parametrized by rational functions of their parameters, they can be viewed as algebraic varieties. This part of OSCAR provides functionality for working with the following commonly used models:

- undirected and directed Gaussian graphical models
- undirected and directed discrete graphical models
- phylogenetic Markov models on trees and networks 
- discrete and Gaussian conditional independence models

Most of the above models are *graphical models* which are parametric statistical models $\mathcal{M}_G$ associated to a graph $G$. In all of the cases described above, the model $\mathcal{M}_G = \phi_{G}$ where $\phi_G$ is a rational map. Key features include:

- Constructing and working with the parametrizations $\phi_G$
- Computing a (partial) generating set of $\text{ker}(\phi_G)$
- Determining the local and global conditional independence statements implied by $G$
- Constructing critical ideals and computing maximum likelihood and euclidean distance degrees



General textbooks offering details on theory and algorithms include:

- [Sul18](@cite)
- [DSS09](@cite)


## Contact

Please direct questions about this part of OSCAR to the following people:
- [Benjamin Hollering](https://sites.google.com/view/benhollering)
- [Marina Garrote López](https://sites.google.com/view/marinagarrotelopez)
- [Tobias Boege](https://taboege.de/)
