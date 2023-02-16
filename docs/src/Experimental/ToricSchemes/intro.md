```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

One can associate to a toric variety a scheme. The resulting class of schemes then benefits from
the functionality that is available both for toric varieties and for schemes.


## Content

We currently provide the following basic conversions:
* Turn an affine toric variety into an affine scheme.
* Turn a normal toric variety into a covered scheme.
Some (but not yet all) of the toric functionality is available for the resulting toric schemes.


## Status

This project is *experimental*. This means that names of functions may change without being
announced. In addition, the  existing functionality has not yet been tested for a long time, and
so may break unexpectedly.


## Long term goals

We aim for a robust interface between toric varieties and schemes.


## Contact

Please direct questions about this part of OSCAR to the following people:
* [Martin Bies](https://martinbies.github.io/),
* [Matthias Zach](https://www.mathematik.uni-kl.de/en/agag/people/members/seite).
