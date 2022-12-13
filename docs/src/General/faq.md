```@contents
Pages = ["faq.md"]
```

# Frequently Asked Questions


## General questions

**Q: How do I install OSCAR?**

You can find our installation instructions [here](https://oscar.computeralgebra.de/install/).

---

**Q: Why do some of your types have funny names like `fmpz` or `fmpq_mat`?**

This has historical reasons. We plan to rename these types before OSCAR 1.0
(the old names will still work indefinitely, though).

---

**Q: Can I find all methods that apply to a given object?**

Yes, Julia provides the function [methodswith](https://docs.julialang.org/en/v1/stdlib/InteractiveUtils/#InteractiveUtils.methodswith) for this very purpose.

For your convenience, let us give an example here. To this end, we first create a projective space in OSCAR:
```julia
julia> v = projective_space(NormalToricVariety,2)
A normal, non-affine, smooth, projective, gorenstein, fano, 2-dimensional toric variety without torusfactor

julia> typeof(v)
NormalToricVariety
```
Suppose that we now want to find all methods that accept a `NormalToricVariety` as one of their arguments.
This can be achieved as follows:
```
julia> methodswith(typeof(v))
[1] intersection_form(v::NormalToricVariety) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/CohomologyClasses/special_attributes.jl:101
[2] mori_cone(v::NormalToricVariety) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/attributes.jl:976
[3] nef_cone(v::NormalToricVariety) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/attributes.jl:953
[4] toric_ideal(ntv::NormalToricVariety) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/attributes.jl:510
[5] volume_form(v::NormalToricVariety) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/CohomologyClasses/special_attributes.jl:50
```
Often it can be beneficial to also include supertypes in the search:
```
julia> methodswith(typeof(v), supertypes = true)
```
As of December 2022, this results in a list of 101 functions.

Note that we can also find the constructors, i.e. functions that return an object of type `NormalToricVariety`.
This is possible with the Julia function [methods](https://docs.julialang.org/en/v1/base/base/#Base.methods):
```julia
julia> methods(typeof(v))
# 5 methods for type constructor:
[1] NormalToricVariety(P::Polyhedron) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/constructors.jl:183
[2] NormalToricVariety(PF::PolyhedralFan) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/constructors.jl:155
[3] NormalToricVariety(rays::Vector{Vector{Int64}}, max_cones::Vector{Vector{Int64}}; non_redundant) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/constructors.jl:131
[4] NormalToricVariety(C::Cone) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/constructors.jl:79
[5] NormalToricVariety(polymakeNTV::Polymake.BigObject) in Oscar at /datadisk/Computer/Mathematics_software/PackagesForJulia/Oscar.jl/src/ToricVarieties/NormalToricVarieties/constructors.jl:8
```

---

**Q: Why do you have your own matrix types, and why do they not support the exact same commands as Julia matrices?**

Unfortunately, Julia's matrices and linear algebra cannot be made to work in
our context due to two independent problems:
  - In empty matrices (0 rows or columns) all that is known is the *type* of
    the matrix entries, however for the complex types used in OSCAR, this
    information is not sufficient to create elements, hence `zero(T)` or
    friends cannot work.
  - Many functions (e.g. `det`) assume that all types used embed into the
    real or complex numbers, in Julia `det(ones(Int, (1,1))) == 1.0`, so the
    fact that this is exactly the integer `1`  is lost. Furthermore, more
    general rings cannot be embedded into the reals at all.
  
---

**Q: Why can `zero(T)` for a type `T` not work?**

At least two reasons:
  - The type depends on data that is not a bit-type.
  - Even if it could, it is not desirable. Typical example: computations in
    ``Z/nZ``, so modular arithmetic. If ``n`` is small, then it is tempting to
    define a type `T` depending on ``n``. We actually did this, and tried to
    use this. It did not work well, for various reasons. E.g.:

    A generic algorithmic pattern for problems over the integers is to solve
    them by solving them modulo ``n`` for many ``n``, e.g. chosen as prime
    numbers, and then to combine them. If the type depends on ``n``, then for
    every prime the code gets compiled, thus negating any advantages from the
    use of modular techniqes.

Of course, one could make the ``n`` an additional parameter to all functions
needing it, but then e.g. addition of matrices would have to be implemented
specifically for this case, negating the advantages of generic
implementations.

In OSCAR, the role of the type is split between the actual Julia type and the
`parent`.

---

**Q: What is a `parent`?**

Almost all element-like objects in OSCAR have a parent, i.e., they belong to
some larger structure. For example algebraic numbers belong to a number field,
modular integers belong to a ring ``Z/nZ``, permutations are elements of
permutation groups and so on. The data common to all such elements is
out-sourced to the parent. For a number field for example, the parent contains
the polynomial used to define the field (plus other information).

Given that a type alone is not large enough to contain the data, the parent is
used. Roughly, outside a function signature, a parent replaces the role of the
type. For example, for a ring element `elm` in OSCAR `zero(parent(elm))` works,
even if `zero(typeof(elm))` may not.

---

**Q: How can I install or access custom GAP packages (e.g. unpublished ones)?**

TODO

---

**Q: Why does my program not terminate?**

Many of the algorithms implemented in OSCAR have a very high complexity. Even
if not calling one of these algorithms directly, you may be using it in the
background. Please read our page on [Complex Algorithms in OSCAR].

---
## Windows specific

**Q: How can I install OSCAR on Windows?**

Please follow [the install instructions on our website](https://oscar.computeralgebra.de/install/).

---

**Q: Why does OSCAR require WSL on Windows?**

Several of the OSCAR corner stones originate from Unix-like operating
systems and have no or only limited native support for Windows.

---

**Q: How can I access Linux files from the Explorer?**

Type `\\wsl$` into the Explorer address bar, then press the Enter key.

---

## Linux specific

**Q: Why can't I install OSCAR using the Julia version installed by my package manager?**

Some Linux distributions unfortunately ship crippled versions of Julia by
default, which prevent OSCAR from working. For example the Debian and Ubuntu
Julia packages are missing some files required by OSCAR. In this case, this can
be resolved by also installing the `libjulia-dev` package.

For this reason, we recommend always using the official Julia binaries
available form the Julia website.

---

**Q: What to do if I get an error similar to ```libstdc++.so.6: version `GLIBCXX_3.4.26'```?**

Sometimes installing or updating OSCAR gives the error ```libstdc++.so.6:
version `GLIBCXX_3.4.26'``` or a similar one.

This typically happens when manually installing Julia using the official Julia
binaries from their website. These bundle their own copy of the C++ standard
library, which can lead to trouble if its version differs from the system's C++
library.

As a workaround, you can rename the copy of the C++ library bundled with Julia,
so that the system copy is used. This can be achieved by executing the
following Julia code:
```julia
  path = Libdl.dlpath("libstdc++")
  mv(path,"$path.bak")
```

If for some reason you need to restore the C++ library bundled with Julia, you
can simply rename it back.

**Q: Why does OSCAR fail to precompile when using it with GNU parallel?**

You get errors like the following when trying to run some script using OSCAR
with GNU parallel:
```
  ERROR: LoadError: InitError: ArgumentError: '.../deps/<something>_jll' exists. `force=true` is required to remove '...' before copying.
```

There was a [bug](https://github.com/JuliaLang/julia/issues/35343) in julia
versions before 1.8 that ignored the parent argument for the `tempname`
function when the `TMPDIR` environment variable is set and GNU parallel by
default sets `TMPDIR` to `/tmp`.

Either upgrade to Julia 1.8 or later, or add `ENV["TMPDIR"]=nothing;` to the
beginning of your julia code (before importing / using Oscar).
