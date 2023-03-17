# Developer guide

## Updating the bibliography

When editing `docs/references.bib` please follow the style of the
existing entries. An easy way to do that is to add your new BibTeX entry,
then run [bibtool](http://www.gerd-neugebauer.de/software/TeX/BibTool/en/)
by invoking it as follows from the root directory of the Oscar.jl repository:

    bibtool docs/references.bib -o docs/references.bib

# To-Dos
* The irrelevant ideal, SR ideal, and ideal of linear relations may need to be modified when the exceptional coordinate "e" is included in the blowup. They are currently set up to work when e is eliminated
* Decide whether to stick with global blowups or use charts
* Consolidate notation about sections and line bundles in the documentation
* The Kodaira type function assumes that the singular locus is given by a single coordinate
* The Kodaira type function only works for codimension 1
* Modify the _ambient_space_from_base function to return ring maps from the base to the ambient space and all other constructors/types to appropriately carry that around, then fix the corresponding bit of code in analyze_fibers
