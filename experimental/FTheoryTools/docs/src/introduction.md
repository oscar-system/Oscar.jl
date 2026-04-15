```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Welcome to FTheoryTools

## Overview

**FTheoryTools** is a computational toolkit within the [OSCAR computer algebra system](https://www.oscar-system.org/), designed to assist researchers in working with F-theory models. It focuses on automating and simplifying calculations involving **singular elliptic fibrations**—key geometric objects in F-theory phenomenology.

While the module is tailored for string theorists, it is equally accessible to mathematicians interested in the rich geometry of singular fibrations, even if they are not familiar with F-theory itself.

This page is meant for *end users* of OSCAR, including students and researchers in mathematics and the natural sciences. No background in string theory or theoretical physics is assumed beyond what is needed to understand the geometry of elliptic fibrations. We encourage interested readers to consult the exposition in [Weigand 2018](@cite Wei18) for more background information.

## Why Use FTheoryTools?

F-theory encodes physical data into the structure of singular elliptic fibrations. In these models:

- The **type and location** of singularities relate directly to gauge groups and matter content.
- **Smooth** fibrations typically yield trivial physics and are therefore less interesting for physical applications.

To analyze the geometry effectively, model builders look for a **crepant resolution** of the singular space—one that preserves the Calabi-Yau condition and retains physical meaning. These resolutions are challenging to compute, especially in higher codimension or for non-toric singularities.

FTheoryTools aims to:

- Automate and streamline the crepant resolution process,
- Make computations more reproducible,
- Extract physically relevant features from the resolved space,
- Offer a modular and extensible framework for researchers in both physics and mathematics.

## Key Features

### Constructing Elliptic Fibrations

FTheoryTools supports construction of elliptic fibrations via:

- [Weierstrass models](@ref),
- [Global Tate models](@ref),
- [Hypersurface models](@ref).

All of these represent singular elliptic fibrations, so many operations and properties are shared across them. This shared functionality is documented at [Functionality for all F-theory models](@ref).

Fibrations can be defined over various base spaces:

- Families of abstract bases,
- Toric varieties (best supported),
- (Planned) General schemes and varieties.

Physically relevant cases often have base dimension 1, 2, or 3, but *FTheoryTools* is not limited to these.

### General Blowups (Beyond Toric)

FTheoryTools enables blowups on arbitrary loci—not just toric centers. This allows users to work with a wider class of singularities, including those without a known toric resolution.

### Literature Models

*FTheoryTools* includes a curated database of well-known F-theory models from the literature. These models are stored using the **MaRDI file format**, a JSON-based format aligned with FAIR data principles:

- **Findability**
- **Accessibility**
- **Interoperability**
- **Reusability**

Learn more about the format:

- [MaRDI project website](https://www.mardi4nfdi.de/about/mission)
- [MaRDI file specification on Zenodo](https://zenodo.org/records/12723387)
- [Design paper](https://link.springer.com/chapter/10.1007/978-3-031-64529-7_25)
- [Serialization docs in OSCAR](https://docs.oscar-system.org/stable/General/serialization/)

Current literature models include:

- [Krause, Mayrhofer, Weigand 2011](@cite KMW12),
- [Morrison, Park 2012](@cite MP12),
- [Lawrie, Schafer-Nameki 2013](@cite LS13),
- [Klevers, Mayorga, Damian, Oehlmann, Piragua, Reuter 2015](@cite KM-POPR15),
- [Cvetič, Klevers, Piragua, Taylor 2015](@cite CKPT15),
- [Taylor, Wang 2015](@cite TW15),
- [Cvetič, Halverson, Ling, Liu, Tian 2019](@cite CHLLT19).

More information: [Literature constructions](@ref).

### ``G_4``-Flux Enumeration

FTheoryTools supports the enumeration of vertical ``G_4``-fluxes—important for understanding chiral spectra in F-theory.

#### Example: [Taylor, Wang 2015](@cite TW15)

- 101 toric rays
- 198 maximal cones
- Defining equation with 355,785 monomials
- 206 toric blowups + 3 smoothness blowups

The computed flux space: ``\mathbb{Z}^{224} \times \mathbb{Q}^{127}``.

Details are available in [BMT25](@cite BMT25). See also: [G4-Fluxes](@ref).

## Tutorials

Explore example-driven tutorials at: [https://www.oscar-system.org/tutorials/FTheoryTools/](https://www.oscar-system.org/tutorials/FTheoryTools/)

These walk through:

- Building models
- Performing blowups
- Extracting physical/geometric information

## Getting Started

1. **Install OSCAR**: [Installation instructions](https://www.oscar-system.org/install/)
2. **Update regularly**: Stay current with new features via the [upgrade guide](https://www.oscar-system.org/upgrade/)

## Project Status

FTheoryTools is an **experimental** module. Most features are well-tested for **toric** models, with active development underway for:

- Support for general base families and schemes
- Line bundle cohomology for vector-like spectra
- Support for terminal/non-minimal singularities
- Base blowup automation

## Contact & Community

For questions, suggestions, or collaboration:

- [Martin Bies](https://martinbies.github.io/)
- [Miķelis E. Miķelsons](https://github.com/emikelsons)
- [Andrew P. Turner](https://apturner.net/)

Community platforms:

- [OSCAR Slack](https://www.oscar-system.org/community/#Slack)
- [GitHub Issue Tracker](https://www.oscar-system.org/community/#Reporting-Issues)

## Acknowledgements

We thank [Mirjam Cvetič](https://live-sas-physics.pantheon.sas.upenn.edu/people/standing-faculty/mirjam-cvetic) and [Mohab Safey El Din](https://www.lip6.fr/actualite/personnes-fiche.php?ident=P816#) for valuable discussions.

The authors are thankful for the support offered by the [TU-Nachwuchsring](https://rptu.de/en/tu-nachwuchsring-network-for-young-scientists-support/home-page). This work was supported by the [_SFB-TRR 195 Symbolic Tools in Mathematics and their Application_](https://www.computeralgebra.de/sfb/) of the [German Research Foundation (DFG)](https://www.dfg.de/en). Martin Bies acknowledges financial support from the \emph{Forschungsinitiative des Landes Rheinland-Pfalz} through the project [_SymbTools -- Symbolic Tools in Mathematics and their Application_](https://fingolfin.github.io/SymbTools/). Andrew P. Turner acknowledges funding from _DOE (HEP) Award DE-SC0013528_ and NSF grant _PHY-2014086_.
