```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
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
linear_program
```

## Solving a linear program - an example
Let $P=[-1,1]^3$ be the $3$-dimensional cube in $\mathbb{R}^3$, and consider
the linear function $\ell$, given by $\ell(x,y,z) = 3x-2y+4z+2$. Minimizing
$\ell$ over $P$ can be done by solving the corresponding linear program.
Computationally, this means first defining a linear program:

```jldoctest
julia> P = cube(3)
Polyhedron in ambient dimension 3

julia> LP = linear_program(P,[3,-2,4];k=2,convention = :min)
Linear program
   min{c⋅x + k | x ∈ P}
where P is a Polyhedron{QQFieldElem} and
   c=Polymake.Rational[3 -2 4]
   k=2
```

The information about the linear program `LP` can be easily extracted.

```jldoctest
julia> P = cube(3);

julia> LP = linear_program(P,[3,-2,4];k=2,convention = :min);

julia> c, k = objective_function(LP)
(QQFieldElem[3, -2, 4], 2)

julia> P == feasible_region(LP)
true
```

To solve the optimization problem call `solve_lp`, which returns a pair `m, v`
where the optimal value is `m`, and that value is attained at `v`.

```jldoctest
julia> P = cube(3);

julia> LP = linear_program(P,[3,-2,4];k=2,convention = :min);

julia> m, v = solve_lp(LP)
(-7, QQFieldElem[-1, 1, -1])

julia> ℓ = objective_function(LP; as = :function);

julia> ℓ(v) == convert(QQFieldElem, m)
true
```

!!! note "Infinite solutions"
    Note that the optimal value may be $\pm\infty$ which currently is
    well-defined by `Polymake.jl`, but not with the `QQFieldElem` number type. Hence
    manual conversion is necessary, until this issue has been resolved.

The optimal value and an optimal vertex may be obtained individually as well.

```jldoctest
julia> P = cube(3);

julia> LP = linear_program(P,[3,-2,4];k=2,convention = :min);

julia> M = optimal_value(LP)
-7

julia> V = optimal_vertex(LP)
3-element PointVector{QQFieldElem}:
 -1
 1
 -1
```


## Functions

```@docs
feasible_region(lp::LinearProgram)
objective_function(lp::LinearProgram{T}; as::Symbol = :pair) where T<:scalar_types
solve_lp(LP::LinearProgram)
optimal_value(lp::LinearProgram{T}) where T<:scalar_types
optimal_vertex(lp::LinearProgram{T}) where T<:scalar_types
load_lp(file::String)
save_lp(target::Union{String,IO}, lp::Union{MixedIntegerLinearProgram{QQFieldElem},LinearProgram{QQFieldElem}})
load_mps(file::String)
save_mps(target::Union{String,IO}, lp::Union{MixedIntegerLinearProgram{QQFieldElem},LinearProgram{QQFieldElem}})
```
