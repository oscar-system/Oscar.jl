# GAP Integration

This section explains how OSCAR interacts with GAP.

## The Julia package [GAP.jl](https://github.com/oscar-system/GAP.jl)

This package provides a bidirectional interface between GAP and Julia.
Its [documentation](https://oscar-system.github.io/GAP.jl/stable/)
describes how to call GAP functions in Julia code and vice versa,
and how low level Julia objects can be converted to GAP objects
and vice versa.

When one works interactively in an OSCAR session,
calling `GAP.prompt()` opens a GAP session which has access to the variables
in the Julia session, in particular to all OSCAR functions and objects;
one can return to the Julia prompt by entering `quit;` in the GAP session.

## Using GAP from OSCAR

Only a part of GAP's functionality has a counterpart in OSCAR.
Everything else can be called directly, which is described here.
Whenever an OSCAR function for the task exists, prefer it: it takes care
of the conversions, and it returns OSCAR objects.

### Global GAP variables and functions

Every global GAP variable is available as `GAP.Globals.<name>`,
and GAP functions can be called like Julia functions.
Julia values of type `Int64` and `Bool` can be passed as arguments
directly; integers of other types, and all other arguments, must be
converted with `GapObj` first, see below.
Small integers, booleans, and elements of small finite fields that GAP
returns are converted to Julia automatically.
Most other results are `GapObj`s and get printed with the prefix `GAP:`
(but a GAP function may also return any Julia object, for example one
that was passed to it as an argument).
```jldoctest
julia> GAP.Globals.Factorial(5)
120

julia> GAP.Globals.IsPrimeInt(97)
true

julia> x = GAP.Globals.Factorial(30)
GAP: 265252859812191058636308480000000
```

### Converting GAP objects to OSCAR objects

For numbers, lists, and strings, use the constructor of the Julia type
that you want to get.
Nested lists need the `recursive` keyword argument, and OSCAR matrices are
constructed from GAP matrices with `matrix`.
```jldoctest
julia> ZZ(GAP.Globals.Factorial(30))
265252859812191058636308480000000

julia> Vector{Int}(GAP.Globals.DivisorsInt(12))
6-element Vector{Int64}:
  1
  2
  3
  4
  6
 12

julia> Vector{Vector{Int}}(GAP.Globals.Partitions(4))
5-element Vector{Vector{Int64}}:
 [1, 1, 1, 1]
 [2, 1, 1]
 [2, 2]
 [3, 1]
 [4]

julia> String(GAP.Globals.StructureDescription(GAP.Globals.SymmetricGroup(4)))
"S4"

julia> matrix(ZZ, GAP.Globals.IdentityMat(2))
[1   0]
[0   1]
```

For example, this is how to call GAP's `Derangements`.
```jldoctest
julia> x = GAP.Globals.Derangements(GapObj(1:4))
GAP: [ [ 2, 1, 4, 3 ], [ 2, 3, 4, 1 ], [ 2, 4, 1, 3 ], [ 3, 1, 4, 2 ], [ 3, 4, 1, 2 ], [ 3, 4, 2, 1 ], [ 4, 1, 2, 3 ], [ 4, 3, 1, 2 ], [ 4, 3, 2, 1 ] ]

julia> Vector{Vector{Int}}(x)
9-element Vector{Vector{Int64}}:
 [2, 1, 4, 3]
 [2, 3, 4, 1]
 [2, 4, 1, 3]
 [3, 1, 4, 2]
 [3, 4, 1, 2]
 [3, 4, 2, 1]
 [4, 1, 2, 3]
 [4, 3, 1, 2]
 [4, 3, 2, 1]
```

GAP groups are turned into OSCAR groups by the constructors `PermGroup`,
`PcGroup`, `FPGroup`, and `matrix_group`.
```jldoctest
julia> G = PermGroup(GAP.Globals.MathieuGroup(11))
Permutation group of degree 11 and order 7920

julia> H = PcGroup(GAP.Globals.SmallGroup(8, 3))
Pc group of order 8
```

Elements of GAP's finite fields and cyclotomic fields are converted by
calling the OSCAR field on them, and the same works for matrices.
```jldoctest
julia> F = GF(3);

julia> F(GAP.Globals.Z(3))
2

julia> K, z = cyclotomic_field(5);

julia> K(GAP.Globals.E(5))
z_5

julia> matrix(K, GAP.evalstr("[[E(5), 1], [0, E(5)^2]]"))
[z_5       1]
[  0   z_5^2]
```

For more complex situations, `Oscar.iso_gap_oscar` and
`Oscar.iso_oscar_gap` construct an isomorphism between a GAP ring or
field and its OSCAR counterpart,
see the section on [GAP Integration](@ref).
More conversions are described in the
[GAP.jl manual](https://oscar-system.github.io/GAP.jl/stable/conversion/).

### Passing OSCAR objects to GAP

The function `GapObj` converts Julia and OSCAR objects into GAP objects.
Not every OSCAR group is based on a GAP group, but permutation groups,
pc groups, finitely presented groups, and matrix groups are;
for those, `GapObj(G)` returns the underlying GAP group.
The same holds for group elements, character tables,
and other objects that are backed by GAP.
```jldoctest
julia> G = symmetric_group(4);

julia> GAP.Globals.StructureDescription(GapObj(G))
GAP: "S4"

julia> GAP.Globals.IsAbelian(GapObj(G))
false

julia> GapObj(ZZ(2)^100)
GAP: 1267650600228229401496703205376

julia> GapObj([1, 2, 3])
GAP: [ 1, 2, 3 ]

julia> GapObj([[1, 2], [3, 4]]; recursive = true)
GAP: [ [ 1, 2 ], [ 3, 4 ] ]

julia> GapObj(matrix(QQ, [1 2; 3 4]))
GAP: [ [ 1, 2 ], [ 3, 4 ] ]
```

Julia functions can be passed to GAP functions that expect a function
argument.
```jldoctest
julia> GAP.Globals.Filtered(GapObj(1:10), is_prime)
GAP: [ 2, 3, 5, 7 ]
```

### Running GAP code

`GAP.evalstr` evaluates a string containing GAP code and returns the
result.
The string has to be a self-contained piece of GAP code;
GAP code that is stored in a file can be read with `GAP.Globals.Read`.
```jldoctest
julia> GAP.evalstr("List([1..5], x -> x^2)")
GAP: [ 1, 4, 9, 16, 25 ]

julia> r = GAP.evalstr("rec(a := 1, b := [1, 2])");

julia> r.a
1
```

`GAP.prompt()` opens a GAP prompt inside the Julia session.
In it, the variables of the Julia session are available via `Julia.<name>`,
and the OSCAR functions via `Oscar_jl.<name>`.
Enter `quit;` to get back to the Julia prompt.
GAP's help system is available at the Julia prompt via
`?GAP.Globals.Size` or `GAP.show_gap_help("Size", true)`.

GAP packages that are not loaded automatically can be loaded with
`GAP.Packages.load`,
and packages that are not shipped with OSCAR can be installed with
`GAP.Packages.install`.
```jldoctest
julia> GAP.Packages.load("ctbllib")
true
```

### When OSCAR does not have the function you need

If you cannot find a counterpart of a GAP function in OSCAR
(the tables in the [Notes for GAP users](@ref) may help),
calling the GAP function as described above will usually work,
as long as the objects involved can be converted.
In this case, please
[open an issue](https://github.com/oscar-system/Oscar.jl/issues)
to request a proper OSCAR interface for it.

Note that `Oscar.GAPWrap`, described below, is intended for OSCAR's own
code, not for interactive use.

## Interface functionalities beyond GAP.jl

For code involving Julia types that are defined in OSCAR,
GAP.jl cannot provide utility functions such as conversions to and from GAP.

- The GAP package OscarInterface (at `gap/OscarInterface`)
  is intended to contain the GAP code in question,
  for example the declarations of new filters
  and the installation of new methods.

  Note that such code must be loaded at runtime into the GAP session
  that is started by Julia, and the OscarInterface package gets loaded
  in OSCAR's `__init__` function.

- The files in the directory `src/GAP`
  are intended to contain the Julia code in question,
  for example conversions from GAP to `ZZRingElem`, `QQFieldElem`,
  `FinFieldElem`, etc.,
  and the construction of isomorphisms between algebraic structures
  such as rings and fields in GAP and OSCAR,
  via [`Oscar.iso_oscar_gap`](@ref) and [`Oscar.iso_gap_oscar`](@ref).

- In OSCAR code, global GAP variables can be accessed as members of
  `GAP.Globals`, but for the case of GAP functions,
  it is more efficient to use `Oscar.GAPWrap` instead.

  For example, if one wants to call GAP's `IsFinite` then it is
  recommended to replace the call `GAP.Globals.IsFinite(x)::Bool`,
  for some GAP object `x` (a group or a ring or a list, etc.),
  by `Oscar.GAPWrap.IsFinite(x)`.
  This works only if the method in question gets defined in
  `src/GAP/wrappers.jl`, thus methods with the required signatures
  should be added to this file when they turn out to be needed.

  (The reason why we collect the `GAP.@wrap` lines in an OSCAR file and
  not inside GAP.jl is that we can extend the list without waiting for
  releases of GAP.jl.)

  Note that `Oscar.GAPWrap` is intended only for *calling* the GAP function
  in question.
  In situations where a GAP function is used for other purposes,
  usually as an argument in a function call,
  one should access it via `GAP.Globals`.

- In GAP code, global Julia variables can be accessed as members of
  `Julia`, relative to its `Main` module.
  For example, one can call `Julia.sqrt` and `Julia.typeof`
  (or `Julia.Base.sqrt` and `Julia.Core.typeof`) in GAP code.

  In order to access variables from the `Oscar` module,
  it is not safe to use `Julia.Oscar`
  because the module `Oscar` is not always defined in `Main`.
  Instead, there is the global GAP variable `Oscar_jl`.

```@docs
Oscar.iso_oscar_gap
Oscar.iso_gap_oscar
```
