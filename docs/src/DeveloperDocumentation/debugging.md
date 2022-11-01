# Debugging OSCAR Code

## Pitfalls: Mutable objects in OSCAR code

Suppose you are having the following difficulties. Your code is exhibiting
inexplicable behaviour and values that should not be changing are changing in
seemingly random locations. To get to the bottom of these kind of issues
it is necessary to be familiar with mutable objects in Julia and some of the
relevant conventions in place in OSCAR. This section discusses these
informal rules as well as some of the exceptions to these rules.

In Julia, objects that can change after construction are declared with the
`mutable struct` keywords and satisfy the `ismutable` predicate.
These objects can be linked together into an arbitrary dependency graph, and
a change to one object may therefore have unintended consequences on another
object in the system.

The simplest example is the creation of a polynomial ring. If we mutate the
array of symbols used for printing, we have effectively changed the ring.

```
julia> v = [:x, :y, :z]; R = PolynomialRing(QQ, v)[1]
Multivariate Polynomial Ring in x, y, z over Rational Field

julia> v[3] = :w; R
Multivariate Polynomial Ring in x, y, w over Rational Field
```

In this example, the modification of `v` is unexpected and may in fact corrupt
the internal data structures used by the polynomial ring. As such, this
modification of `v` has to be considered illegal. Upon creation of the array
called `v`, we have full rights over the object and can mutate at will.
However, after passing it to the function `PolynomialRing`, we have given up
*owernership* of the array and are no longer free to modify it.

General OSCAR Principle (GOP):

*Code should be expected to behave as if all objects are immutable.*

Ramifications:

1. This means that the polynomial ring constructor is allowed to expect that `v`
   is never mutated for the remaining duration of its life. In return, the
   constructor is guaranteed not to modify the array, so that `v` is still
   `[:x, :y, :z]` after `PolynomialRing` returns.
2. In general this means that all functions should be expected to take ownership
   of their arguments: the user is safest *never* modifying an existing object
   that has been passed to an unknown Julia function. Note that assignments
   such as `a[i] = b` or `a.foo = b` usually mutate the object `a`.
   See [Ownership of function arguments](@ref)
3. For reasons of efficiency, it is sometimes desirable to defy this principle
   and modify an existing object. The fact that a given function may modify a
   preexisting object is usually communicated via coding conventions on the
   name - either a `!` or a `_unsafe` in the name of the function.
   See [Unsafe arithmetic with OSCAR objects](@ref)

### Ownership of function arguments

In this example we construct the factored element `x = 2^3` and then change the
`2` to a `1`. The GOP says this modification of `a` on line 3 is illegal.

```julia
julia> a = fmpz(2)
2

julia> x = FacElem([a], [fmpz(3)]); evaluate(x)
8

julia> a = one!(a)  # illegal in-place assignment of a to 1
1

julia> evaluate(x)  # x has been changed and possibly corrupted
1
```

In the previous example, the link between the object `x` and the object `a` can
be broken by passing a `deepcopy` of `a` to the `FacElem` function.

```julia
julia> a = fmpz(2)
2

julia> x = FacElem([deepcopy(a)], [fmpz(3)]); evaluate(x)
8

julia> a = one!(a)  # we still own a, so modification is legal
1

julia> evaluate(x)  # x is now unchanged
8
```

It is of course not true that all Julia functions take ownership of their
arguments, but the GOP derives from the fact that this decision is an
implementation detail with performance consequences. The behaviour of a
function may be inconsistent across different types and versions of OSCAR.
In the following two snippets, the GOP says both modifications of `a` are
illegal since they have since been passed to a function. If `K = QQ`, the two
mutations turn out to be legal currently, while they are illegal if
`K = quadratic_field(-1)[1]`. Only with special knowledge of the types can the
GOP be safely ignored.

```
R = PolynomialRing(K, [:x, :y])[1]
a = one(K)
p = R([a], [[0,0]])
@show p
a = add!(a, a, a)       # legal? (does a += a in-place)
@show p
```

```
R = PolynomialRing(K, :x)[1]
a = [one(K), one(K)]
p = R(a)
@show (p, degree(p))
a[2] = zero(K)          # legal?
@show (p, degree(p))
```

### Ownership of function return values

The nuances of who is allowed to modify an object returned by a function is
best left to the next section [Unsafe arithmetic with OSCAR objects](@ref).
The GOP says of course you should not do it, but there are cases where it can
be more efficient. However, there is another completely different issue of
return values that can arise in certain interfaces.

First, we create the Gaussian rationals and the two primes above `5`.

```julia
julia> K, i = quadratic_field(-1)
(Imaginary quadratic field defined by x^2 + 1, sqrt(-1))

julia> m = Hecke.modular_init(K, 5)
modular environment for p=5, using 2 ideals
```

The function `modular_project` returns the projection of an element of `K` into
each of the residue fields.

```julia
julia> a = Hecke.modular_proj(1+2*i, m)
2-element Vector{fq_nmod}:
 2
 0
```

While the function has produced the correct answer, if we run it again on a
different input, we will find that `a` has changed.

```julia
julia> b = Hecke.modular_proj(2+3*i, m)
2-element Vector{fq_nmod}:
 1
 3

julia> a
2-element Vector{fq_nmod}:
 1
 3
```

The preceeding behaviour of the function `modular_proj` is an artifact of
internal efficiency and may be desirable in certain circumstances. In other
circumstances, the following `deepcopy`s may be necessary for your code to
function correctly.


```julia
julia> a = deepcopy(Hecke.modular_proj(1+2*i, m));
julia> b = deepcopy(Hecke.modular_proj(2+3*i, m));
julia> (a, b)
(fq_nmod[2, 0], fq_nmod[1, 3])
```

## Unsafe arithmetic with OSCAR objects

Particularly with integers (`BigInt` and `fmpz`) - but also to a lesser extent
with polynomials - the cost of basic arithmetic operations can easily be
dominated by the cost of allocating space for the answer. For this reason,
OSCAR offers an interface for in-place arithmetic operations.

Instead of writing `x = a + b` to compute a sum, one writes `x = add!(x, a, b)`
with the idea that the object to which `x` is pointing is modified instead of
having `x` point to a newly allocated object. In order for this to work, `x`
must point to a *fully independent* object, that is, an object whose
modification through the interface [Unsafe operators](@ref) will not change
the values of other existing objects.
The actual definition of "fully independent" is left to the implementation of
the ring element type. For example, there is no distinction for immutables.

It is generally not safe to mutate the return of a function. However, the
basic arithmetic operations `+`, `-`, `*`, and `^` are guaranteed to return
a fully independent object regardless of the status of their inputs.
As such, the following implementation of `^` is illegal by this guarantee.

```julia
function ^(a::RingElem, n::Int)
  if n == 1
    return a    # must be return deepcopy(a)
  else
    ...
  end
end
```

In general, if you are not sure if your object is fully independent,
a `deepcopy` should always do the job.
