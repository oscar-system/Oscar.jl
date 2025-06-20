```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Welcome to FTheoryTools

## Goal

*FTheoryTools* is a computational toolkit designed to automate and simplify many of the intricate and time-consuming computations encountered in F-theory model building. Developed as a submodule of the `OSCAR` computer algebra system, it provides infrastructure to streamline model construction, crepant resolution of singularities, and geometric analysis of singular elliptic fibrations—central objects in F-theory phenomenology.

F-theory models encode physical data in the geometry of singular elliptic fibrations. The presence and type of singularities in the fibration are physically meaningful, determining the resulting gauge symmetry and matter content. Smooth fibrations typically correspond to trivial physics and are not the focus of model builders. However, performing computations directly on singular spaces is often impractical. Therefore, F-theory studies usually begin by finding a *crepant resolution* of the singular geometry—one that preserves the Calabi-Yau condition and hence the physical relevance.

The crepant resolution of singularities, especially in higher codimensions or non-toric settings, poses major computational and conceptual challenges:
1. Not all singularities admit crepant resolutions.
2. When crepant resolutions exist, they may not resolve all singularities.
3. Multiple resolutions might exist, offering complementary physical insights.

While *FTheoryTools* does not (yet) resolve all such foundational issues, its aim is to *make the crepant resolution process as accessible, automated, and reproducible as possible*. It also assists in extracting physically relevant geometric features from the resolved space.


## Key Features

### Elliptic Fibration Models

You can construct elliptic fibrations using:
- [Weierstrass models](@ref),
- [Global Tate models](@ref),
- [Hypersurface models](@ref).
Since all of these different models encode singular elliptic fibrations mathematically, certain properties, attributes, and methods apply to all models equally. This shared functionality is summarized on the page [Functionality for all F-theory models](@ref).

The base of the fibration may be:
- A family of abstract base spaces,
- A toric variety (the currently best-supported case),
- In the future: schemes and more general varieties.

We do not limit ourselves to base spaces of dimensions 1, 2, or 3, but those are of particular interest for physical applications.

### General Blowups

Unlike most traditional approaches which are limited to toric blowups, *FTheoryTools* supports blowups on arbitrary loci, including non-toric ones. This significantly expands the scope of F-theory constructions accessible through our software.

### Literature Models

We provide support for an expanding database of well-known F-theory models, referred to as *literature models*. These models are serialized using the **MaRDI file format**, which is based on JSON and designed in accordance with the FAIR principles --- Findability, Accessibility, Interoperability, and Reusability.

The MaRDI file format is part of a broader effort by the [Mathematical Research Data Initiative (MaRDI)](https://www.mardi4nfdi.de/about/mission), a German consortium focused on establishing standards and developing tools to improve the handling of mathematical research data. While OSCAR's current implementation of the MaRDI format is centered around saving and loading data within the system, the format is designed to be extensible and interoperable with other computer algebra systems such as `SAGE` and `Macaulay2`.

The specification of the **mrdi** file format used in OSCAR can be found on [Zenodo](https://zenodo.org/records/12723387). For a broader discussion of its design and goals, refer to [this article](https://link.springer.com/chapter/10.1007/978-3-031-64529-7_25), as well as the [OSCAR documentation on serialization](https://docs.oscar-system.org/stable/General/serialization/).

Current database entries include, but are not limited to, models from:
- [Krause, Mayrhofer, Weigand 2011](@cite KMW12),
- [Morrison, Park 2012](@cite MP12),
- [Lawrie, Schafer-Nameki 2013](@cite LS13),
- [Klevers, Mayorga, Damian, Oehlmann, Piragua, Reuter 2015](@cite KM-POPR15),
- [Cvetič, Klevers, Piragua, Taylor 2015](@cite CKPT15),
- [Taylor, Wang 2015](@cite TW15),
- [CHLLT 2019](@cite CHLLT19).

More details are available on the page [Literature constructions](@ref).

### ``G_4``-Flux Enumeration

*FTheoryTools* supports enumeration of ``G_4``-fluxes, essential for determining the chiral spectra of F-theory compactifications. As a demonstration of its capabilities, we did enumerate ``G_4``-fluxes for the notoriously difficult *F-theory geometry with the most flux vacua* [Taylor, Wang 2015](@cite TW15), involving:
- 101 toric rays,
- 198 maximal cones,
- A defining equation with 355,785 monomials,
- A resolution via 206 toric blowups (plus 3 additional ones for ambient smoothness).
We identified a space of candidate vertical ``G_4`` fluxes parametrized by ``\mathbb{Z}^{224} \times \mathbb{Q}^{127}``. Details on this particular computation are available in the publication [BMT25](@cite BMT25).

More details are available on the page [G4-Fluxes](@ref).

## Tutorials

We encourage you to take a look at the tutorials on FTheoryTools, which can be found at [https://www.oscar-system.org/tutorials/FTheoryTools/](https://www.oscar-system.org/tutorials/FTheoryTools/).


## Installation

To get started with *FTheoryTools*, install [OSCAR](https://www.oscar-system.org/install/).

To benefit from the latest features, improvements, and bug fixes, we recommend keeping OSCAR up to date by following the [upgrade instructions](https://www.oscar-system.org/upgrade/).


## Status and Outlook

*FTheoryTools* is an experimental module. Functionality is currently most robust for model based on toric geometries, but support for general families and schemes is actively being developed. Future extensions will include, but are not limited to:
- Vector-like spectrum computations via line bundle cohomologies,
- Support for terminal and non-minimal singularities,
- Automated base blowups.


## Contact

For questions or feedback, reach out to:
- [Martin Bies](https://martinbies.github.io/)
- [Miķelis E. Miķelsons](https://github.com/emikelsons)
- [Andrew P. Turner](https://apturner.net/)

Or visit:
- [OSCAR Slack](https://www.oscar-system.org/community/#Slack)
- [GitHub Issue Tracker](https://www.oscar-system.org/community/#Reporting-Issues)


## Acknowledgements

We thank [Mirjam Cvetič](https://live-sas-physics.pantheon.sas.upenn.edu/people/standing-faculty/mirjam-cvetic) and [Mohab Safey El Din](https://www.lip6.fr/actualite/personnes-fiche.php?ident=P816#) for valuable discussions.

The authors are thankful for the support offered by the [TU-Nachwuchsring](https://rptu.de/en/tu-nachwuchsring-network-for-young-scientists-support/home-page). This work was supported by the [_SFB-TRR 195 Symbolic Tools in Mathematics and their Application_](https://www.computeralgebra.de/sfb/) of the [German Research Foundation (DFG)](https://www.dfg.de/en). Martin Bies acknowledges financial support from the \emph{Forschungsinitiative des Landes Rheinland-Pfalz} through the project [_SymbTools -- Symbolic Tools in Mathematics and their Application_](https://fingolfin.github.io/SymbTools/). Andrew P. Turner acknowledges funding from _DOE (HEP) Award DE-SC0013528_ and NSF grant _PHY-2014086_.
