```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
```

# Mixed Integer Linear Programs


## Introduction

The purpose of a mixed integer linear program is to optimize a linear function
over a polyhedron, where we require a given subset of the variables to be
integers.



## Constructions

Mixed integer linear programs are constructed from a polyhedron and a linear
*objective function* which is described by a vector and (optionally) a
translation. One can select whether the optimization problem is to maximize or
to minimize the objective function. Furthermore one can optionally specify
`integer_variables`, a subset of variables that are required to be integral,
given as a set of integers, their respective indices.

```@docs
mixed_integer_linear_program
```

## Functions
```@docs
feasible_region(milp::MixedIntegerLinearProgram)
ambient_dim(milp::MixedIntegerLinearProgram)
optimal_value(milp::MixedIntegerLinearProgram{T}) where T<:scalar_types
optimal_solution
solve_milp
integer_hull
gomory_chvatal_closure
```
