```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Welcome to FTheoryTools

## Goal

We aim to automate numerous recurring and, at least in part, tedious computations in F-Theory model building, as
detailed extensively in [Wei18](@cite). The primary focus of our software is on (complex) elliptic fibrations with
singularities. Those singularities hold significant importance. In essence, the absence of these singularities
implies that the geometry encodes trivial or uninteresting physics. Therefore, the smooth case is typically not explored.

A substantial amount of information about the physics is encoded in the geometry of this singular fibration. To some
extent it is clear what geometric quantities are to be considered, to some extent this is an open question. Regardless,
computing quantities of interest, such as intersection theory and Chern classes, on singular spaces is challenging. However,
it is possible to link the physics encoded by the singular geometry to the physics encoded by related smooth geometries.
Consequently, almost all F-Theory studies begin by computing a smooting-out, i.e. a resolution, of the singular geometry in
question.

In order to easily link the physics on the singular geometry to the physics of one of its resolution, the resolution in question
must be crepant, i.e. must preserve the Calabi-Yau condition. However, this restriction to crepant resolutions introduces additional challenges:
1. A crepant resolution cannot resolve all singularities, meaning certain singularities may persist. This is an area of interest in recent F-Theory investigations.
2. The existence of a crepant resolution is not guaranteed.
3. Exploring whether different resolutions provide insights into different aspects of the singular geometry poses further questions. However, just as it is unclear whether a single crepant resolution is known, there is currently no way to confirm that all crepant resolutions have been identified.

*FTheoryTools* may not (yet) answer these profound and fundamental questions. Instead, the goal of this suite of computer tools is to streamline and simplify the crepant singularity resolution process as much as possible, as well as the subsequent extraction of geometric
features of the resolved space.

With *FTheoryTools*, you can create elliptic fibrations using one of the following models:
* Weierstrass model,
* Global Tate model,
* Hypersurface model.
In each case, the base space can be a family of spaces (as used in the literature when an explicit base is not specified) or a toric space. We anticipate an extension to schemes as bases. The dimension of those bases is not limited to $1$, $2$, or $3$, but of course includes those cases important to the physics.

We also offer a database of models frequently studied in the literature. For those models we use the term 'literature_models.' In essence, you should be able to create a model described in a paper at the click of a button. We anticipate that this feature will greatly simplify future research in F-Theory.


## Status

We anticipate the following workflow:

### User Input:
- Create your desired F-Theory model by using one of the methods described above.
- Choose a resolved phase/crepant resolution.
- Select generating sections for $\operatorname{U}(1)$ symmetries.

### Output:
- (Crepantly) resolved geometry.
- Singular loci in suitable codimension (e.g., $1$, $2$, and $3$ if the base space has dimension $3$).
- Fiber diagrams of the resolved fiber over the originally singular loci, including intersections of $\operatorname{U}(1)$-sections.
- Gauge group.
- Topological data (e.g., Euler number).

Currently, our primary focus is on the Elephant in the room, that is the crepant resolution. While already functional for numerous setups, especially literature models, it is far from complete. At this point, the largest functionality exists for toric cases, with ongoing work to extend those toric resolution techniques to families of spaces and schemes.

We are also actively expanding our database of supported literature models. At this point, amongst others, our database includes models from the following papers:
- The Tate Form on Steroids: Resolution and Higher Codimension Fibers [LS13](@cite),
- F-Theory on all Toric Hypersurface Fibrations and its Higgs Branches [KM-POPR15](@cite),
- Quadrillion F-Theory Compactifications with the Exact Chiral Spectrum of the Standard Model [CHLLT19](@cite).

Consequently, we are already covering infinite families of models (e.g., with $SU(k)$ gauge group, where $k \geq 1$).
Still, the total number of papers in our database is at this point limited to about $6$. This is to be extended a lot.


## Tutorial

We encourage you to take a look at the tutorials on FTheoryTools, which can be found
[here](https://www.oscar-system.org/tutorials/FTheoryTools/).


## Possible future extensions

Future extensions include, but are not necessarily limited to, the following:
* Specify a $G_4$-flux and work out the chiral spectra,
* Specify a gauge potential and work out (candidates for) the line bundles whose cohomologies encode the vector-like spectra,
* Other singularity types (non-minimal, terminal, etc.,)
* Base blowups for singularity resolution.


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Martin Bies](https://martinbies.github.io/),
* [Mikelis Emils Mikelsons](https://github.com/emikelsons),
* [Andrew Turner](https://apturner.net/).

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).
Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).


## Acknowledgements

We appreciate insightful discussions with [Mirjam Cvetič](https://live-sas-physics.pantheon.sas.upenn.edu/people/standing-faculty/mirjam-cvetic) and 
[Mohab Safey El Din](https://www.lip6.fr/actualite/personnes-fiche.php?ident=P816#). Martin Bies and Mikelis Mikelsons appreciate support by the TU-Nachwuchsring. The work of Andrew Turner is supported by DOE (HEP) Award DE-SC001352.
