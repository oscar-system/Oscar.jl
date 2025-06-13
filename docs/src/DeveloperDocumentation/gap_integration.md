# GAP Integration

This section explains how Oscar interacts with GAP.

## The Julia package [GAP.jl](https://github.com/oscar-system/GAP.jl)

This package provides a bidirectional interface between GAP and Julia.
Its [documentation](https://oscar-system.github.io/GAP.jl/stable/)
describes how to call GAP functions in Julia code and vice versa,
and how low level Julia objects can be converted to GAP objects
and vice versa.

When one works interactively in an Oscar session,
calling `GAP.prompt()` opens a GAP session which has access to the variables
in the Julia session, in particular to all Oscar functions and objects;
one can return to the Julia prompt by entering `quit;` in the GAP session.

## Interface functionalities beyond GAP.jl

For code involving Julia types that are defined in Oscar,
GAP.jl cannot provide utility functions such as conversions to and from GAP.

- The GAP package OscarInterface (at `gap/OscarInterface`)
  is intended to contain the GAP code in question,
  for example the declarations of new filters
  and the installation of new methods.

  Note that such code must be loaded at runtime into the GAP session
  that is started by Julia, and the OscarInterface package gets loaded
  in Oscar's `__init__` function.

- The files in the directory `src/GAP`
  are intended to contain the Julia code in question,
  for example conversions from GAP to `ZZRingElem`, `QQFieldElem`,
  `FinFieldElem`, etc.,
  and the construction of isomorphisms between algebraic structures
  such as rings and fields in GAP and Oscar,
  via [`Oscar.iso_oscar_gap`](@ref) and [`Oscar.iso_gap_oscar`](@ref).

- In Oscar code, global GAP variables can be accessed as members of
  `GAP.Globals`, but for the case of GAP functions,
  it is more efficient to use `Oscar.GAPWrap` instead.

  For example, if one wants to call GAP's `IsFinite` then it is
  recommended to replace the call `GAP.Globals.IsFinite(x)::Bool`,
  for some GAP object `x` (a group or a ring or a list, etc.),
  by `Oscar.GAPWrap.IsFinite(x)`.
  This works only if the method in question gets defined in
  `src/GAP/wrappers.jl`, thus methods with the required signatures
  should be added to this file when they turn out to be needed.

  (The reason why we collect the `GAP.@wrap` lines in an Oscar file and
  not inside GAP.jl is that we can extend the list without waiting for
  releases of GAP.jl.)

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
