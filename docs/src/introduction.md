```@meta
CurrentModule = Oscar
```

# Gröbner walk

This is a sample text to outline the structure of the packages in the `Experimental` folder.
You can show docstrings like this:

The general idea of the Gröbner walk was proposed by [AGK97](@cite).
```@docs
    groebner_walk(
      I::MPolyIdeal, 
      target::MonomialOrdering = lex(base_ring(I)),
      start::MonomialOrdering = default_ordering(base_ring(I));
      perturbation_degree = ngens(base_ring(I)),
      algorithm::Symbol = :standard
    )
```
    
