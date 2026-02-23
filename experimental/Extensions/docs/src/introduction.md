```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Extensions

An extension is a module that is loaded when all of its dependencies are loaded in the current Julia session.
This functionality has been available since Julia version 1.9. 
See the [julia docs](https://docs.julialang.org/en/v1/manual/code-loading/#man-extensions) on how to setup an extension.

We will maintain extensions of Oscar so long as the dependencies used as part of the extension are actively maintained

## Homotopy Continuation

Homotopy continuation is a general concept that comprises of methods for numerically solving systems of polynomial equations.
[HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) is a Julia package for working with such methods.

Our extension allows users to call functions from HomotopyContinuation.jl on Oscar types.
This is done by overloading HomotopyContinuation.jl contructors `Expression`, `System` and `witness_set`. 

To align with Oscar naming conventions and to raise awareness that the underlying methods are numerical we have exported the following functions.

```@docs
solve_numerical
dim_numerical
```

## Contact

Please direct questions about this part of OSCAR to the following people:
* [Antony Della Vecchia](https://antonydellavecchia.github.io/)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
