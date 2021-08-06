```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_linear_programs.md"]
```

# Linear programs


## Introduction

The task of a linear program is to optimize a linear function over a polyhedron.



## Constructions

Linear programs are constructed by passing a polyhedron, and a linear function in terms of a vector and translation. One may specify whether the optimization problem is in terms of maximum or minimum.

```@docs
LinearProgram
```

## Solving a linear program

After constructing a linear program, it is straightforward to extract the feasible region and objective function

```@docs
feasible_region
```

```@docs
objective_function
```

Calling `solve_lp` on a linear program outputs a pair: the optimal value and the vertex at which the optimum is obtained

```@docs
solve_lp
```

One can obtain the maximal or minimal value of the objective function over the feasible region, if it exists.

```@docs
maximal_value
```
```@docs
minimal_value
```

One can also obtain the maximal or minimal vertex at which the objective function attains its maximal or minimal value (respectively).

```@docs
maximal_vertex
```

```@docs
minimal_vertex
```
