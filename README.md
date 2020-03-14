# Oscar.jl

Welcome to the OSCAR project, a visionary new computer algebra system
which combines the capabilities of four cornerstone systems: Gap,
Polymake, Antic and Singular.

## Installation

* Download Julia version 1.3 or higher
* In Julia type
  - using Pkg
  - Pkg.add("Oscar")
  - using Oscar

Installation will take about a minute.

## Current requirements

* we support Julia 1.3 and we try to get ready a.s.a.p. for the upcoming Julia 1.4
* target platforms: x86_64-linux-gnu, x86_64-apple-darwin14
* CxxWrap 0.9.0
* travis: distribution ubuntu bionic

Some explanations:

* Julia 1.0 is the latest LTS version, but it seems unlikely that GAP (with its garbage collection) can ever be supported
* Julia 1.1 and 1.2 do not receive any back ports any more
* Windows support only through WSL (i.e., this is covered by x86_64-linux-gnu)
* The travis configuration is important for building CxxWrap; otherwise the compiler is too old.


## Examples of usage

.... coming soon.
