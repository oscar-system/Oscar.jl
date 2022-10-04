# Design Decisions
This document covers the ideas and design decisions behind Oscar, as well as
some pitfalls to avoid.


## Oscar - what is the idea
Oscar is the innovative, next generation Computer Algebra System. The ultimate
goal for Oscar is to compete with (and ideally beat) Magma and Sage in our
areas of expertise. Oscar should be accessible, even for the youngest student
who is familiar with these objects. Oscar should follow general mathematical
conventions to support the widest possible range of applications.

The key idea for development of OSCAR is to pick a single textbook for every
area and follow the conventions in there.  This way we get a consistent
interface.

## Oscar and Julia

Oscar is written in Julia, but is not Julia, nor can it be. 
Some examples to illustrate what that means:
Julia's matrices are arrays (of arbitrary dimension), parameterized by
the type of the entries (apart from banded, sparse, ... special matrices).
In the numerical world, the type mostly defines the representation of an object
- double and variations
- complex
- BigFloat
- Int
- BigInt

In algebra, this is either not true or terribly inefficient (or impossible)
Take $\mathbb{Z}/n\mathbb{Z}$ integers modulo $n$, and matrices over it
   
Then either:
-  $n$ is part of the type -> every function is recompiled for every $n$ -
  which kills all modular (Chinese remainder theorem  (CRT) based)
  algorithms
-  $n$ is not part of the type, then it needs to be elsewhere, e.g. in the
   parent, or in every element, or by passing additional arguments, or ...

Furthermore, if $n$ is BigInt (fmpz), so no bittype, then it cannot be part of
the type.

For non-empty matrices, this can be compensated if the entries store enough
information, but for empty matrices this information needs to be collected
elsewhere.

To summarize, normal Julia infrastructure does not suffice for our purposes in
many places. Hence we provide our own, which any code contributions should use.
If functions are missing in it, then please
- add them
- or tell us

## Mathematical Context in Oscar

When studying mathematics, the exact meaning of a term or object is determined
by context. In Oscar, this does not work. The meaning has to be part of either

 - the object
 - or the question posed about the object
 
As an example: in classical number theory there is the convention that many
definitions that are trivial for fields are silently applied to the ring of
integers. One speaks of the unit group of the number field, meaning the unit
group of the ring of integers. Let alpha be an element explicitly constructed
as an element of the number field, not of the ring of integers, then

- `is_unit` will just test if it is non-zero (unit in a field as a special type
  of ring)
- `is_unit_in_ring_of_integers` would supply the context for the other
  interpretation.

In Oscar, this context is mostly supplied by the type of the object and
possibly the parent, e.g. the containing ring/ field/ group.

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
For specialized functionality we can access the entirety of the following
software frameworks on a lower level:
 - polymake
 - Singular
 - Gap

So: Please use it. It is safe to assume all can be improved, however, if we try
to perfect every single line of code again and again, we we won't get anywhere;
there is a balance to be found.
For preference: 
```
Correctness > Interoperability, Readability
Interoperability > Speed
Readability > Speed
```
Having said that: of course, sometimes pure speed matters, but not nearly as
often as people think.


## What Is Missing?


The infrastructure is incomplete, e.g. we do not have
 - varieties, (we've started with curves recently)
 - combinatorial manifolds
 - surfaces
 - tropical polytopes
 - ....

Since we are a relatively small team and OSCAR is still very new, the usual

    I work for 6 month in a separate repo on a branch and then will dazzle you
    with perfect code and cool examples

approach is not going to work for now. It will result in everyone fixing the same
infrastructure problems over and over again. Please consider to work, e.g. in
a file/directory  in `Oscar/examples` and push on a regular basis, even, or in
particular, incomplete code. Break it down into small pull requests. Please also see
the [Introduction for new developers](@ref).



## Practical Development

### Creating New Basic Types
If you encounter the need for a new basic type, say a new multivariate
ring, please consider the ramifications:
 - can you do matrices?
 - modules?
 - ideals?
 - graded stuff?
 - "complete" arithmetic?
 - interaction with other types? (map to residue rings, apply automorphisms, ...)
 - in fact, everything the other MPoly type can?

If no: at least use "our" types to interface your function, better still, use
our type and complain about lack of functionality/ speed/ interface (or provide
patches).

This applies to all foundations! They are all incomplete, and they can all
be improved BUT if everyone does their own foundations, we cannot work
together.

!!! danger "Expert definitions"
    As a reminder, please stick to "global definitions" and not "experts"
    versions of definitions. Reasoning such as: "but all experts know and
    expect this - it is always done this way" will make your function
    impossible to be used by outsiders. Feel free to add the other "expert"
    definition layer if you need.

### Adding new dependencies
Unless really really necessary, don't add new dependencies. Every new
dependency complicates the development workflow, in that we will need to stay
compatible with this package. 

