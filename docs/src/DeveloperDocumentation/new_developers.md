# Introduction for new developers

This document is meant to get new developers started. It will not go into depth
of programming in Julia or working with git, as there are far better resources
on these things online.

!!! note "Pay attention to your GitHub notifications!"
    Once you open a pull request on GitHub you will receive feedback, comments,
    and questions on GitHub. So please pay attention to your GitHub
    notifications.

## Oscar - what is the idea

Oscar is the innovative, next generation Computer Algebra System. The ultimate goal for Oscar is to compete with (and ideally beat) Magma and Sage in our areas of expertise. Oscar should be accessable, even for the youngest student who is familiar with these objects. Oscar should follow general mathematical conventions to support the widest possible range of applications.

## Oscar and Julia

Oscar is written in Julia, but is not Julia, nor can it be. 
 Some examples to illustrate what that means:
  Julias matrices are (C) arrays (of arbitrary dimension), parametrized by
  the type of the entries (apart from banded, sparse, ... special matrices).
  In the numerical world, the type mostly defines the representation of an object
  - double and variations
  - complex
  - BigFloat
  - Int
  - BigInt

  In algebra, this is either not true or terribly inefficient (or impossible)
  Take
    Z/nZ integers modulo n, and matrices over it
   
   Then either
    n is part of the type -> every function is recompiled for every n - which kills
    all modular (Chinese remainder theorem  (CRT) based) algorithms

    n is not part of the type, then it need to be elsewhere, ie. in the parent,
    every element, additional arguments, ...

    Furthermore, if n is BigInt (fmpz), so no bittype, then it cannot be part of the type.

  for non-empty matrices, this does not matter as the entries have enough information. For
  empty matrices....

  So: offshot: we simply cannot use normal Julia infrastructure in many places, so please
  use ours. If functions are missing then
  - add them
  - or tell us

## Mathematical Context in Oscar

When studying mathematics, the exact meaning of a term or object is determined by context. In 
Oscar, this does not work: the meaning has to be part of either

 - the object
 - or the question posed about the object
 
As an example: in classical number theory there is the convention that many definitions
that are trivial for fields are silently applied to the ring of integers. One
speaks of the unit group of the number field, meaning the unit group of the ring
of integers. Let alpha be an element explicitly constructed as an element of the number field, not of the ring of integers, then

   - is_unit will just test if it is non-zero (unit in a field as a special type of ring)
   - is_unit_in_ring_of_integers would supply the context for the other interpretation.

In Oscar, this context is mostly supplied by the type of the object and possibly the 
parent, e.g. the ring/ field/ group we're in.

## What Do We Have:

We have a large codebase for infrastructure in place, comprising at least
 - matrices
 - polynomials (univariate and multivariate)
 - power series
 - number fields
 - (abelian) groups
 - polytopes, cones, linear programs
 - polyhedral fans
 - ... and MUCH more
As well as specialised stuff
 - polymake
 - Singular
 - Gap

So: please use it. It is safe to assume all can be improved, however, if we try to perfect every
single line of code again and again, we we won't get anywhere; there is a balance to be found.
For preference: interoperability is much more
important than speed, readability is much more important than speed. Having
said that: of course, sometimes pure speed matters, but not nearly as often as
people think.

The infrastructure is incomplete, e.g. we do not have
 - varieties, (we've started with curves recently)
 - simplicial complexes, combinatorial manifolds
 - surfaces
 - tropical varieties, in particular hypersurfaces
 - tropical polytopes
 - ....

If you need them, e.g. you want a library/ package for (elliptic) curves, then
step ONE should be to define a curve. It does not need to be perfect, support
all of Hartshorn immediately, but it also should not contradict him either. It
is OK for a function to say: sorry, not (yet) implemented. Then if I need it,
I can add and re-use the infrastructure.

The key idea is to pick a single textbook and follow the conventions in there. This
way we get a consistent interface.

## What Is Missing?

At this point we have, mathematically speaking, almost nothing. This has to
change, quickly.  We work in a team of ~ 10(?) people.
We lack infrastructure as everyone will notice immediately.
As a result, the classical
   I work for 6 month in a separate repo on a branch and then will dazzle you with perfect
   code and cool examples
is not going to work for now. It will result in everyone fixing the same
infrastructure problems over and over again. Please consider to work, e.g. in
a file/directory  in `Oscar/examples` and push on a regular basis, even, or in
particular, incomplete code. Break it down into small pull requests.
 - it is in examples, so no-one is accidentally going to automatically use it
 - it is in the main repo, so everyone can see it and try to use it, producing feedback
 - infrastructure, in Oscar at least, can be added in the same repo, immediately and thus
   is available to others, immediately. (feedback, ...)
 - Talk to us - as early as possible. There is a chance that what you're trying to
   do is already there - in a different package maybe, not linked in, ....
 - Talk to us!
   - on slack: send an email to webmaster-oscar@mathematik.uni-kl.de to be added
   - on oscar-dev@mathematik.uni-kl.de
   - whenever you meet us (at the uni, during conferences, ...)

## Important notes
1. If you encounter error messages after rebasing to the current master chances
   are that some dependencies need upgrading. Please first try whether
   executing `]up` gets rid of your errors.
2. Please have a look at the [Developer Style Guide](@ref). Adhering to the
   style guide makes reviewing code easier for us, and hence your new feature
   can be merged faster.
3. Let us know what you are working on:
   - You can open a draft pull request on GitHub right at the beginning of your
     work.
   - Feel free to contact us on
     [Slack](https://join.slack.com/t/oscar-system/shared_invite/zt-thtcv97k-2678bKQ~RpR~5gZszDcISw).
   - Have a look at [our community page](https://oscar.computeralgebra.de/community/).
4. Please also read our page on [Documenting OSCAR code](@ref).
5. Look at existing code that does similar things to your project to get an
   idea of what OSCAR code should look like. Try to look at multiple examples.

## Overview
In general you have to do
the following six steps for submitting changes to the Oscar source:

1. Fork the main Oscar.jl repository. For this go to the [Oscar.jl GitHub
   page](https://github.com/oscar-system/Oscar.jl) and click on "Fork" at the
   upper right.
2. Clone your forked repository to your local machine.
   ```
   git clone git@github.com:your_github_username/Oscar.jl
   ```
3. Create a new branch, usually the naming convention is to use your initials
   ("yi") and then describe your change, for example:
   ```
   git checkout -b yi/new_feature
   git checkout -b yi/issue1234
   git checkout -b yi/document_feature
   ```
4. Edit your source and try out your changes locally (see below). To use your local copy of
   the sources, start Julia and
   ```
   ]dev /path/to/local/clone/of/your/fork/of/Oscar.jl
   ```
   If this succeeds, you can enter `using Oscar` in Julia and it will use your local
   copy.
5. Once you are done editing, push your branch and open a pull request. It is
   recommended that you open a draft pull request to the main Oscar repository
   as soon as you start working. That way Oscar developers are aware of work
   being done and can give feedback early in the process.
6. Once you have finished your work, mark your pull request as ready. It will
   then be reviewed and, probably after feedback and requests for changes,
   merged.

## Alternative: `]dev Oscar`

Alternatively you can call
```
]dev Oscar
```
in Julia. This will create a directory `~/.julia/dev/Oscar`. This directory is
a git clone of the central Oscar repository. You can develope your code here,
however you will still have to fork Oscar, as you have no rights to push to the
central repository. You can then add your fork as another remote, have a look
at the section on rebasing below for hints.

## Practical Development
### Developing Efficiently
For practical development:
 - use [the `Revise` package](https://github.com/timholy/Revise.jl)
 - start a file/ directory in Oscar/examples, call it Blah.jl

Then `Oscar.example("Blah.jl")` will include it, so you don't have to worry about a path;
and `Oscar.revise("Blah.jl")` will put it under revise, i.e. re-load it automatically whenever you make a change.

If Blah.jl starts with

```julia
 module BlahModule
 using Oscar

 ...

 end

 using .BlahModule
``` 

then you can even redefine struct's without quitting Julia (and restarting)

### Creating New Basic Types
Whenever you feel the urge to define a new basic type, say a new multivariate
ring, please consider the ramifications:
 - can you do matrices?
 - modules?
 - ideals?
 - graded stuff?
 - "complete" arithmetic?
 - interaction with other types? (map to residue rings, apply automorphisms, ...)
 - in fact, everything the other MPoly type can?

if yes: still reconsider, but go ahead

if no: at least use "our" types to interface your function, better still, use our type and
complain about lack of functionality/ speed/ interface (or provide patches)

This applies to all foundations! They are all incomplete, and they can all
be improved BUT if everyone does their own foundations, we cannot work
together.

As a reminder, please stick to "global definitions" and not "experts" versions of definitions. Reasoning such as: "but all experts know and expect this - it is always done this way" will make your function impossible to be used by outsiders. Feel free to add the other "expert" definition layer if you
need.

### Including Outside Packages
Unless really really necessary, don't include outside packages
as dependencies. Every outside package means that the testing, releases, ...
are more complicated.  You're at the whim of others.  If you need a highly
non-trivial 10000 line package, well, there is no help and we need to see how
to arrange it. In particular, for the near future, it will be your job to
ensure that this outside package still works with us.  If you need a trivial 10
line function from a 10000 line package, just don't.

## The edit process
### Editing the source
The sources can be found in the `src` folder. Please pay attention to the
folder structure and choose sensibly where to place your code (when fixing a
bug this is probably a minor question).

### Adding tests

### Adding documentation
There are two places where documentation can be added:
1. In the docstrings above the functions in the `src` folder;
2. In the documentation files in the `docs/src` folder.
3. The overall structure is fixed in the file `docs/doc.main`. If you create a
   new file in `docs/src`, you will have to add an entry in `docs/doc.main`.

In general, 1 is preferred to 2, i.e. any explanation of the functions and
objects should go there and the files in `docs/src` should remain relatively
sparse. Please also pay attention to the documentation section of the
[Developer Style Guide](@ref).


## Further hints

### Ask Oscar related questions in the Oscar slack

### Use `]up`
Working with the development version also entails that the packages Oscar
depends on need to be up to date. Julia can update these packages if you type
`]up` in the Julia prompt. Many error messages after updating the source can be
resolved by simply updating.

### Style guide
Please have a look at the [Developer Style Guide](@ref) to get an overview over
naming conventions, code formatting, etc.

### Building the documentation
To build and test the documentation, please have a look at [Documenting OSCAR
code](@ref).

### Rebasing
One way to stay up to date with the current master is rebasing. In order to do
this, add the main Oscar.jl repository as a remote, fetch, and then rebase.
```
git remote add oscar-system git@github.com:oscar-system/Oscar.jl
git fetch oscar-system
git rebase oscar-system/master
```
Adding the remote only has to be executed once.

### Pitfalls: Mutable objects in Oscar code

Suppose you are having difficulties with your code: it is exhibiting
inexplicable behaviour and values that should not be changing are changing in
seemingly random locations. To get to the bottom of these kind of issues
it is necessary to be familiar with mutable objects in Julia and some of the
relevant conventions in place in Oscar. This section discusses these
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

General Oscar Principle (GOP):

*Code should be expected to behave as if all mutable objects are immutable.*

Ramifications:
1. This means that the polynomial ring constructor is allowed to expect that `v`
   is never mutated for the remaining duration of its life. In return, the
   constructor is guaranteed not to modify the array, so that `v` is still
   `[:x, :y, :z]` after the constructor returns.
2. In general this means that all functions should be expected to take ownership
   of their arguments: the user is safest *never* modifying an existing object
   that has been passed to an unknown Julia function. Note that assignments
   such as `a[i] = b` or `a.foo = b` usually mutate the object `a`.
   See [Ownership of function arguments](@ref)
3. For reasons of efficiency, it is sometimes desirable to defy this principle
   and modify an existing object. The fact that a given function may modify a
   preexisting object is usually communicated via coding conventions on the
   name - either a `!` or a `_unsafe` in the name of the function.
   See [Unsafe arithmetic with Oscar objects](@ref)

#### Ownership of function arguments

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
function may be inconsistent across different types and versions of Oscar.
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

#### Ownership of function return values

The nuances of who is allowed to modify an object returned by a function is
best left to the next section [Unsafe arithmetic with Oscar objects](@ref).
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

#### Unsafe arithmetic with Oscar objects

Particularly with integers (`BigInt` and `fmpz`) - but also to a lesser extent
with polynomials - the cost of basic arithmetic operations can easily be
dominated by the cost of allocating space for the answer. For this reason,
Oscar offers an interface for in-place arithmetic operations.

Instead of writing `x = a + b` to compute a sum, one writes `x = add!(x, a, b)`
with the idea that the object to which `x` is pointing is modified instead of
having `x` point to a newly allocated object. In order for this to work, `x`
must point to a *fully independent* object, that is, an object whose
modification through the interface [Unsafe operators](@ref) will not change
the values of other existing objects.
The actual definition of "fully independent" is left to the implementation of
the ring element type. For example, there is no distinction for immutables.

It is generally not safe to mutate the return value of a function. However, the
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
