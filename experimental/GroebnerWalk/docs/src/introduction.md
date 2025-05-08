```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Usage

The Gröbner walk is an approach to reduce the computational complexity of Gröbner basis computations as proposed by [AGK97](@cite).
These incarnations of the Gröbner walk refer to a family of algorithms that perform a reverse local search on the cones of the Gröbner fan.
Then, a Gröbner basis is calculated for each encountered cone while reusing the generators obtained from the previous cone.

The implemented algorithms may be accessed using the following function.

```@docs
    groebner_walk(
      I::MPolyIdeal, 
      target::MonomialOrdering = lex(base_ring(I)),
      start::MonomialOrdering = default_ordering(base_ring(I));
      algorithm::Symbol = :standard
    )
```
    
