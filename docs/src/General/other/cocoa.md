# Notes for CoCoA users

These notes are aimed at people who have experience using CoCoA-5, and
are considering trying out OSCAR.  As always, the passage from one
system that one knows well to a new system can be challenging -- we hope
that the notes here will be helpful.

Importantly, in tandem with these notes, the reader should also consult
the section about Julia's foibles, especially the handling of literal
values: see [Integers and rational numbers](@ref other_integers).

The guidance here should also benefit CoCoALib users contemplating
experimenting with OSCAR, but we do not give hints about the passage
from C++ to Julia.


## Why switch to OSCAR?

OSCAR offers a much wider range of built-in functions than CoCoA, though
not every function in CoCoA has a counterpart in OSCAR.  In many cases,
the implementations in OSCAR are more advanced, so they may handle a wider
range of inputs, and may be significantly faster on larger inputs,
though there are also some cases where CoCoA is comparably fast.

It makes sense to consider using OSCAR if you have reached CoCoA's limits.
You could also contact the maintainers of CoCoA to see whether they can
help alleviate the limits -- but you might need to be patient!


## How to obtain/install OSCAR

The main OSCAR website can be found at https://www.oscar-system.org/
The menu on that page offers a link to `Installation` which explains
the steps necessary to download and install OSCAR -- you will need a
good internet connection, and some patience to complete the installation.
Be aware that the installation requires as much as 8Gbyte disk space;
heed also the warnings at https://www.oscar-system.org/install/linux/


## Forewarning: just-in-time compilation and garbage collection

_Forewarned is forearmed:_  your first impression of OSCAR/julia could
easily be that it is terribly slow!  This is largely illusory.

An important difference between CoCoA and OSCAR is that CoCoA is
pre-compiled C++ code (so starts up quickly), whereas OSCAR runs in a
Julia environment which performs so-called _just-in-time compilation_
(aka. _stop-and-make-the-user-wait_).  This manifests itself in two ways:
initial start-up of OSCAR may take several seconds; and the first time
that an OSCAR function is called, it will be compiled by Julia, which
may produce a noticeable "hesitation".  However, once the function has
been compiled, it will likely run faster next time it is called.

Another consequence of just-in-time compilation is that the Julia
session quickly grows to occupy 2Gbyte or more RAM: this is largely
workspace for the Julia compiler and the compiled OSCAR code.  After
this initial surge in memory consumption, further demand for RAM is
more moderate.

Julia uses a garbage-collecting memory manager.  This means that Julia
will occasionally stop running the program for a short while so that
the garbage collector can do its work.  This can also manifest itself
as brief "hesitations".


## Naming conventions in OSCAR

CoCoA-5 uses mostly _CamelCase_ for the names of functions, and
some longer words are abbreviated: for example, `IsPolyRing`.

In contrast OSCAR uses _snake_case_ for the names of functions,
and abbreviations are rarely used: for example, `is_polynomial_ring`.


## Julia language and syntax

OSCAR is written in Julia, and to use OSCAR you must use the Julia
language.  We cannot put a Julia language tutorial here, as it would
be either too large or woefully incomplete.  We mention a few _gotchas_
of Julia: _e.g._ commands whose meaning or effect may be unexpected.


### Semicolons

Commands in CoCoA-5 end with a semicolon (or an `End...` keyword).
So a command in CoCoA may be spread over several lines.

In Julia a command may end with a semicolon or with a newline.  In
general, the semicolon is used to _suppress printing_ of the result
in interactive use.  A long command can be spread over several lines
thanks to the following rule:  if the command is syntactically
incomplete when newline is reached then the next line is interpreted
as part of the command.


### Assignment and equality testing

In CoCoA-5, to assign a value to a variable one uses `:=`; in Julia
one uses just `=`.

In CoCoA-5, the operator `=` tests equality, while the operator `<>`
tests for non-equality.  In Julia, the equality test uses `==` (double
equals sign), and non-equality uses `!=`.

Julia's conventions are like their counterparts in C++.


### Integers

In CoCoA-5 all integers are of unlimited precision.  In OSCAR/Julia the
true situation is more complicated, and can easily catch out the unwary!
Be sure to read carefully the explanation in
[Integers and rational numbers](@ref other_integers) 
In summary: integer literals work fine as indices, for other uses it
is safer to write the literal inside a function call `ZZ(...)`.


### Division and rationals

We cannot change it, but Julia's interpretation of division of integers
is most inconvenient in a computer algebra context.

Actually, Julia offers two division operators, namely `/` and `//`.
These are explained in [/ vs //](@ref subtle_distinction_for_rings).

In summary: **Don't write a rational as `1/3`**; the safest alternative
is to write `QQ(1,3)`.  As a compromise, `1//3` might suffice in many
circumstances.


### Current ring

In CoCoA-5 there is the notion of _current ring_ to facilitate input
of polynomials in an interactive session.  OSCAR does not have a
corresponding mechanism.  Nevertheless, unless you frequently change
the current ring in CoCoA-5, the OSCAR approach will work fine for you.

In both CoCoA-5 and OSCAR you have to create a polynomial ring before
you can compute with its elements -- see the examples below, and look
at the OSCAR documentation for its polynomial rings:
https://docs.oscar-system.org/stable/AbstractAlgebra/mpolynomial/#polynomial_ring-Tuple%7BRing,%20Int64%7D

OSCAR offers a much wider range of rings than CoCoA-5 does.


### Using output as input

In CoCoA-5, most printed values can be copied and used as valid input.
In OSCAR this works in only a few cases (but the OSCAR printed form
is often neater than that produced by CoCoA-5).


## Some examples

### Factorizing a polynomial

Observe that in both cases we must create the polynomial ring before
giving the polynomial as input, and then operating on it.

In CoCoA-5:
```
/**/ use P ::= ZZ[x];
/**/ f := 4*x^8 + 4*x^6 + x^4 + 4*x^2 + 4;
/**/ indent(factor(f));
record[
  RemainingFactor := 1,
  factors := [2*x^4 +4*x^3 +5*x^2 +4*x +2,  2*x^4 -4*x^3 +5*x^2 -4*x +2],
  multiplicities := [1,  1]
]
```

The same computation in OSCAR:
```jldoctest
julia> P, x = ZZ["x"];

julia> f = 4*x^8 + 4*x^6 + x^4 + 4*x^2 + 4;

julia> factor(f)
(2*x^4 - 4*x^3 + 5*x^2 - 4*x + 2) * (2*x^4 + 4*x^3 + 5*x^2 + 4*x + 2)


### Matrices

The Julia system has its own notion of matrices, but it is simplest to avoid
these when doing computer algebra computations.  This example illustrates
how Fibonacci numbers can be computed via matrix powers; the main purpose
is to show how to create a matrix -- the built-in `fibonacci` function is
a much better way to compute these numbers!

Both CoCoA-5 and OSCAR come with the rings `ZZ` and `QQ` predefined, so
we may use them immediately; moreover, their names are the same in the
two systems.

In CoCoA-5:
```
/**/  M := mat(ZZ, [[0,1], [1,1]]);
/**/  M^25;
matrix(ZZ,
 [[46368, 75025],
  [75025, 121393]])
```

In OSCAR:
```jldoctest
julia> M = matrix(ZZ, [[0,1], [1,1]]);

julia> M^25   # no semicolon means print the result
[46368    75025]
[75025   121393]
```
Observe that the columns are neatly aligned; the lack of commas is a
mixed blessing.


### Groebner bases

We recall that a Groebner basis is a set of polynomial ideal generators
with some nice properties, and the definition of Groebner basis depends
on a term ordering (aka. monomial ordering).  In CoCoA-5 the term ordering
is specified when the polynomial ring is created (defaulting to `degrevlex`).
In OSCAR, the term ordering is specified when the function `groebner_basis`
is called (again defaulting to `degrevlex`).  Here is an example to
illustrate the similarities and differences between CoCoA-5 and OSCAR.

In CoCoA-5 we compute `degrevlex` and `lex` Groebner bases of "the same" ideal:
```
/**/ use P ::= QQ[x,y,z]; // default term order degrevlex
/**/ I := ideal([x^2-y, x*y-z]);
/**/ ReducedGBasis(I);
[x*y -z,  x^2 -y,  y^2 -x*z]

/**/ use PP ::= QQ[x,y,z],lex;
/**/ J := ideal([x^2-y, x*y-z]);
/**/ ReducedGBasis(J);
[x*y -z,  x^2 -y,  x*z -y^2,  y^3 -z^2]

```
We can see that the ideals `I` and `J` are obviously the same, but for
CoCoA-5 they are different because they belong to different polynomial
rings (which are obviously canonically isomorphic).


In OSCAR we can do the same computations a bit more simply:
```jldoctest
julia> R, (x, y, z) = polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal(R, [x^2 - y, x*y - z]);

julia> groebner_basis(I)
Gröbner basis with elements
  1: y^2 - x*z
  2: x*y - z
  3: x^2 - y
with respect to the ordering
  degrevlex([x, y, z])

julia> groebner_basis(I; ordering=lex(R))
Gröbner basis with elements
  1: y^3 - z^2
  2: x*z - y^2
  3: x*y - z
  4: x^2 - y
with respect to the ordering
  lex([x, y, z])
```
In OSCAR we needed to create just one polynomial ring, and just the single
ideal `I`.

