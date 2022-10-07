```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["linear_programs.md"]
```

# Linear Programs


## Introduction

The purpose of a linear program is to optimize a linear function over a
polyhedron.



## Constructions

Linear programs are constructed from a polyhedron and a linear *objective function*
which is described by a vector and (optionally) a translation. One can select whether the optimization
problem is to maximize or to minimize the objective function.

```@docs
LinearProgram
```

## Solving a linear program - an example
Let $P=[-1,1]^3$ be the $3$-dimensional cube in $\mathbb{R}^3$, and consider
the linear function $\ell$, given by $\ell(x,y,z) = 3x-2y+4z+2$. Minimizing
$\ell$ over $P$ can be done by solving the corresponding linear program.
Computationally, this means first defining a linear program:

```@repl oscar
P = cube(3)
LP = LinearProgram(P,[3,-2,4];k=2,convention = :min)
```

The information about the linear program `LP` can be easily extracted.

```@repl oscar
c, k = objective_function(LP)
P == feasible_region(LP)
ℓ = objective_function(LP; as = :function)
```

To solve the optimization problem call `solve_lp`, which returns a pair `m, v`
where the optimal value is `m`, and that value is attained at `v`.

```@repl oscar
m, v = solve_lp(LP)
ℓ(v) == m
```

The optimal value and an optimal vertex may be obtained individually as well.

```@repl oscar
M = optimal_value(LP)
V = optimal_vertex(LP)
ℓ(V) == M
```


## Functions

```@docs
feasible_region(lp::LinearProgram)
objective_function(lp::LinearProgram{T}; as::Symbol = :pair) where T<:scalar_types
solve_lp(LP::LinearProgram)
optimal_value(lp::LinearProgram{T}) where T<:scalar_types
optimal_vertex(lp::LinearProgram{T}) where T<:scalar_types
```

## Saving and loading

Objects of type `LinearProgram` can be saved to a file and loaded from a file
in the following way:
```@repl oscar
C = cube(3)
LP=LinearProgram(C, [1,2,-3], convention=:min)
save("lp.poly", LP)
LP0 = load("lp.poly")
solve_lp(LP0)
solve_lp(LP)
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and Oscar.

