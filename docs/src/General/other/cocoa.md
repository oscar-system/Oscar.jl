# Notes for CoCoA-5 and CoCoALib users

These notes are aimed at people who have experience using CoCoA-5 (and/or
perhaps C++ programs using CoCoALib), and are considering trying out OSCAR.
As always, the passage from one system that one knows well to a new system
can be challenging -- we hope that the notes here will be helpful.


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
Be aware that the installation requires as much as 8Gbyte disk space.
It is also recommended that you use a computer which typically has at
least 8Gbyte RAM free under normal circumstances -- ideally the computer
should have _at least_ 16Gbyte RAM.


## Forewarning: just-in-time compliation

_Forewarned is forearmed:_  your first impression of OSCAR/julia could
easily be that it is terribly slow!  This is somewhat illusory.

An important difference between CoCoA and OSCAR is that CoCoA is
pre-compiled C++ code (so starts up quickly), whereas OSCAR runs in a
Julia environment which performs so-called _just-in-time compilation_
(aka. _stop-and-make-the-user-wait_).  The first time that an OSCAR
function is called, it will be compiled by Julia, and this may produce
a noticeable "hesitation" before the computation begins.  However,
once the function has been compiled, it will run faster next time you
call it.

Another consequence of just-in-time compilation is that the Julia
session quickly grows to occupy 2Gbyte or more RAM: this is largely
workspace for the Julia compiler and the compiled OSCAR code.  After
this initial surge in memory consumption, further demand for RAM is
more moderate.


## Julia language and syntax

OSCAR is written in Julia, and to use OSCAR you must use the Julia
language.  We cannot put a Julia language tutorial here, as it would
be either too large or woefully incomplete.  We mention a few _gotchas_
of Julia: _e.g._ commands whose meaning or effect may be unexpected.


### Semicolons

Commands in CoCoA-5 end with a semicolon (or an `End...` keyword).
So a command in CoCoA may be spread over several lines.

In Julia a command may end with a semicolon or with a newline; in OSCAR
it is customary to write just one command on each line, and use the
newline to demarcate its end.  However, if the command is syntactically
incomplete when newline is reached then the next line is interpreted
as part of the command.


### Assignment and equality testing

In CoCoA-5, to assign a value to a variable one uses `:=`; in Julia
one uses just `=`.  Oddly, if you use `:=` by accident then Julia
issues an absurd error message.

In CoCoA-5 one tests equality using the operator `=`, and non-equality
using the operator `<>`.  In Julia, the equality test uses `==` (double
equals sign), and non-equality uses `!=`.

Julia's conventions are like their counterparts in C++.


### Integers

In CoCoA-5 all integers are of unlimited precision.  In OSCAR/Julia the
true situation is complicated.  Simplifying somewhat, an integer literal
(a number written out as decimal digits) has 64-bit limited precision;
we must write the literal inside a function call `ZZ(...)` to obtain
an unlimited value.  The limited precision integers work well as indexes
into arrays _etc._

OSCAR/Julia does not automatically switch to unlimited integers:
integer overflow is reported only sometimes.  For instance `factorial(99)`
will trigger an error in OSCAR, but `factorial(ZZ(99))` will return
the desired result.  In contrast `99^99` returns a negative result
because it is computed as a signed 64-bit value; `ZZ(99)^99` will
produce the correct answer -- why is `ZZ(99^99)` also wrong?

The CoCoA-5 approach is easier to use, but can sometimes be slow.


### Division and rationals

We cannot change it, but Julia's interpretation of division of integers
is most inconvenient in a computer algebra context.

Actually, Julia offers two division operators, namely `/` and `//`.

**Don't write a rational as `1/3`** because Julia will silently convert
the result to an approximate floating-point value, which you almost
certainly do not want!  In contrast `1//3` does sort of work, but there
are better ways (see below).

OSCAR behaves much better.  Let `a` and `b` be two values.  _If at least
one of then is an OSCAR value_ then `a/b` means what you'd expect: it is
exact division.

If `a` and `b` are both integer values, and at least one of them is
an OSCAR `ZZ` value then `a//b` produces the expected rational number.
Alternatively, you can use `QQ(a,b)` -- this alternative form is likely
safer (and it is obviously a rational number).  Later you can learn
further uses for `//` in OSCAR.


### Current ring

In CoCoA-5 there is the notion of _current ring_ to facilitate input
of polynomials in an interactive session.  OSCAR does not have a
corresponding mechanism.  Nevertheless, unless you frequently change
the current ring in CoCoA-5, the OSCAR approach will work fine for you.

In both CoCoA-5 and OSCAR you have to create a polynomial ring before
can compute with it elements.  OSCAR offers a much wider range of rings
than CoCoA-5 does.


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
```
julia> P,x = ZZ["x"];
julia> f = 4*x^8 + 4*x^6 + x^4 + 4*x^2 + 4;
julia> factor(f)
(2*x^4 - 4*x^3 + 5*x^2 - 4*x + 2) * (2*x^4 + 4*x^3 + 5*x^2 + 4*x + 2)
```


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
Observe that the columns are neatly aligned; the lack of commas is a
mixed blessing.


### Groebner bases

**INCOMPLETE**
