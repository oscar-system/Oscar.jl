# Documenting OSCAR code

The general philosophy of the OSCAR documentation is to put as much of the
information as possible into the docstrings and only use the doc pages for
collecting this information and provide some additional general context.
Exceptions to this philosophy are the developer and general pages.

## Docstrings of exported functions

Exported function should have docstrings, which look like
```julia
Markdown.@doc doc"""
    functionname(x::ArgumentType, b::OtherArgument; c::Keyword = default) -> Int, Int

A short description of the function. It is allowed to use $\LaTeX$.
"""
functionname(x...,b...; c = ...)
```
If the signature is too long, use linebreaks to fit 80 characters.

Please also do provide an example within the docstring if possible, preferrably
as a `jldoctest`, i.e.
````julia
Markdown.@doc doc"""
    functionname(x::ArgumentType, b::OtherArgument; c::Keyword = default) -> Int, Int

A short description of the function. It is allowed to use $\LaTeX$.

# Examples
This shows that `functionname` does the right thing for input `input`
```jldoctest
julia> input = ...

julia> functionname(input)
output
```
"""
functionname(x...,b...; c = ...)
````
This allows the user to immediately see how the function can be used, gives
them some code that they can copy-paste and manipulate, and, as a bonus,
provides a testcase as well.


## The folder `docs`

The folder `docs/src` contains the OSCAR documentation website. Most of the
pages are relatively sparse and consist of
````
```@docs
some_function
some_other_function
[...]
```
````
blocks that simply pull in the docstring from the corresponding source file. If
you add a new page in `docs/src`, you will have to modify `docs/doc.main` to
include your new page in the appropriate place.


## Building the OSCAR documentation

To build the OSCAR documentation (especially when editing it, which
usually requires rebuilding it several times to test the outcome), we
recommend the following steps:

1. Install the Julia package `Revise` into the environment in which you
   run your Oscar dev version:

        using Pkg ; Pkg.add("Revise")

2. Start a fresh Julia session and load Revise before Oscar:

        using Revise, Oscar

3. Build the manual as follows:

        Oscar.build_doc()

4. To rebuild the documentation, just repeat the same command as in step 3,
   without exiting Julia. Thanks to the Revise package, any runs after
   the first will be much faster.

### Running the doctests
For running the doctests, `build_doc` takes an optional parameter `doctest`.
The default value is `false`, so to run the doctests run
```
build_doc(doctest = true)
```
The documentation in the github actions is build using Julia 1.6. Differences
in output may be caused by a difference in Julia versions.

You can also use this parameter to autofix your documentation, in case you
changed the output of some function appearing in many places. See the section
on repairing `jldoctest`s below.


## Updating the bibliography

When editing `docs/oscar_references.bib` please follow the style of the
existing entries. An easy way to do that is to add your new BibTeX entry,
then run [bibtool](http://www.gerd-neugebauer.de/software/TeX/BibTool/en/)
by invoking it as follows from the root directory of the Oscar.jl repository:

    bibtool docs/oscar_references.bib -o docs/oscar_references.bib


## Automatically repairing `jldoctest`s

It is possible to have julia fix the output of all `jldoctest`s when your
changes to the code entail changes to the output. Just run the following
command:
```
build_doc(doctest = :fix)
```
!!! danger
    Please use this command carefully:
    - Make sure to only commit the changes to the doctests originating from
      your changes to the code.
    - The doctests also serve as actual tests, so make absolutely sure that the
      output is still mathematically correct.

!!! tip
    If this command fails with an error message indicating lacking permissions
    to change `AbtractAlgebra.jl` related docs, it may help to run the
    following command:
    ```
    ]dev AbstractAlgebra
    ```
