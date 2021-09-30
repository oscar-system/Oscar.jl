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

Despite the linear program being initially defined with the `:min` convention,
one can still ask for the maximal value and a point on the cube for which
$\ell$ is maximized.

```@repl oscar
M = maximal_value(LP)
V = maximal_vertex(LP)
ℓ(V) == M
```


## Functions

After constructing a linear program, it is straightforward to extract the
feasible region and objective function

```@docs
feasible_region
objective_function
```

Calling `solve_lp` on a linear program outputs a pair: the optimal value and
the vertex at which the optimum is obtained

```@docs
solve_lp
```

One can obtain the maximal or minimal value of the objective function over the
feasible region, if it exists.

```@docs
maximal_value
minimal_value
```

One can also obtain the maximal or minimal vertex at which the objective
function attains its maximal or minimal value (respectively).

```@docs
maximal_vertex
minimal_vertex
```

## Saving and loading

Objects of type `LinearProgram` can be saved to a file and loaded from a file
in the following way:
```@repl oscar
C = cube(3)
LP=LinearProgram(C, [1,2,-3], convention=:min)
save_linearprogram(LP, "lp.poly")
LP0 = load_linearprogram("lp.poly")
solve_lp(LP0)
solve_lp(LP)
```
The file is in JSON format and contains all previously gathered data belonging
to the underlying polymake object. In particular, this file can now be read by
both polymake and Oscar.

```@docs
save_linearprogram(LP::LinearProgram, filename::String)
load_linearprogram(filename::String)
```
