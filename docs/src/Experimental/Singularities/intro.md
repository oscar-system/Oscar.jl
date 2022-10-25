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

The singularities part of OSCAR provides functionality for handling

- space germs
- map germs

For readers' convenience, the documentation is organised by typical settings, 
This, on the other hand, implies that certain keywords may appear in several 
settings. In particular, there are overlaps between hypersurface 
singularities and curve singularities and between hypersurface singularities
and isolated complete intersection singularities

!!! note
    Most of the functions discussed here rely on standard bases. These are implemented in OSCAR for localizations of multivariate polynomial rings over fields (exact fields supported by OSCAR) at points with coordinates in said field. This explained in more detail in [Generalities on Space Germs](@ref space_germ_generalities).


Textbooks offering details on theory (and some algorithms) include:
- [GLS07] (@cite)   
- [dJP00] (@cite)
- [CLS20] (@cite)
- [CLS21] (@cite)

    


