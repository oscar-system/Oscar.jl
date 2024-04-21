# Gröbner walk

GroebnerWalk.jl is a Julia package providing implementations of Gröbner walk algorithms
for computing Gröbner bases over fields on top of Oscar.jl.

## Usage

GroebnerWalk.jl provides its entire functionality through the function `groebner_walk`.
The following example demonstrates the usage. First, we define the ideal Oscar.jl.
```julia
using Oscar

R, (x,y) = QQ[:x, :y]                  # define ring ...
I = ideal([y^4+ x^3-x^2+x,x^4])        # ... and ideal
```

Then, we can pass the ideal to `groebner_walk` to calculate the Gröbner basis.
Here, we want a Gröbner basis with respect to the lexicographic ordering on `R`.
```julia
using GroebnerWalk

groebner_walk(I, lex(R)) # compute the Groebner basis
```

## Status
This repository represents the status of the code as a submission for MEGA 2024.
At the moment, the standard walk (TODO: reference) and the generic walk (TODO: references)
are implemented.
It is slated for inclusion into OSCAR as experimental package.

## Contacts
The library is maintained by Kamillo Ferry (kafe (at) kafe (dot) dev) and Francesco Nowell (francesconowell (at) gmail (dot) com).

## Acknowledgement
The current implementation is based on an implementation by Jordi Welp. We thank him for 
laying the groundwork for this package.
