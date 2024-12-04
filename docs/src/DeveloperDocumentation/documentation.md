# Documenting OSCAR code

The general philosophy of the OSCAR documentation is to put as much of the
information as possible into the docstrings and only use the doc pages for
collecting this information and provide some additional general context.
Exceptions to this philosophy are the developer and general pages.


## Docstrings of exported functions

Exported function should have docstrings, which look like
```julia
@doc raw"""
    functionname(x::ArgumentType, b::OtherArgument; c::Keyword = default) -> Int, Int

A short description of the function. It is allowed to use $\LaTeX$.
"""
functionname(x...,b...; c = ...)
```
If the signature is too long, use linebreaks to fit 80 characters.

Please also do provide an example within the docstring if possible, preferably
as a `jldoctest`, i.e.
````julia
@doc raw"""
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


## Building the OSCAR documentation with `Oscar.build_doc`

!!! note "Previewing the documentation"
    Once you have created a pull request it is possible to preview the
    documentation on github using the link
    https://docs.oscar-system.org/previews/PR<prnumber>/
    where you insert the number of your PR for `prnumber`. Alternatively you
    can look at the github actions tab of your PR and click the details link
    next to the `documenter/deploy` action. There are a few conditions for this 
    to work:
    - No conflicts with the master branch.
    - Documentation action is successful, i.e. no doctest errors.
    - The branch for the PR is in the main `oscar-system/Oscar.jl` repository.
    You can still build the documentation locally with the commands described below.

```@docs
build_doc
```
Please also read the section below on repairing the `jldoctest`s using
`build_doc`.
!!! note "Browser reports denied access"
    Depending on your system, it might happen that the browser opens after a
    successful build, but only informs you that the access to the file was denied.
    This happens, for example, on Ubuntu which comes with a sandboxed Firefox.
    In this case, using `build_doc` with `start_server = true` should circumvent
    this problem.


### Automatically repairing `jldoctest`s

It is possible to have julia fix the output of all `jldoctest`s when your
changes to the code entail changes to the output. Just run the following
command:
```julia
build_doc(doctest = :fix)
```
If you just want to fix some of the `jldoctest`s, and do not want to build
the documentation, you can also use `Oscar.doctest_fix`:
```@docs
Oscar.doctest_fix
```
!!! danger
    Please use these commands carefully:
    - Make sure to only commit the changes to the doctests originating from
      your changes to the code.
    - The doctests also serve as actual tests, so make absolutely sure that the
      output is still mathematically correct.

!!! tip
    If these commands fail with an error message indicating lacking permissions
    to change `AbstractAlgebra.jl` related docs, it may help to run the
    following command:
    ```
    ]dev AbstractAlgebra
    ```


## Updating the bibliography

When editing `docs/oscar_references.bib` please follow the style of the
existing entries. An easy way to do that is to add your new BibTeX entry,
then run [bibtool](http://www.gerd-neugebauer.de/software/TeX/BibTool/en/)
by invoking it as follows from the root directory of the Oscar.jl repository:

    bibtool docs/oscar_references.bib -o docs/oscar_references.bib

For every pull request on github, the CI checks if running `bibtool` leads to
changes in the bibliography. If so, this test fails and indicates that the
(recently) added bibliography entries are not standardized. For a merge, it
is not required that this test is passed. Still, please feel encouraged to fix
this failure by running `bibtool` locally as explained above.

!!! note "bibtool produces changes in unrelated parts of oscar_references.bib"
    Sometimes `bibtool` will produce many changes when run locally. This can be
    caused by a version difference. The version used in our github actions is
    2.68. Check your version by running `bibtool -V`. When running this
    command, please also pay attention whether any "Special configuration
    options" are set.

Please follow the additional guidelines below, that are not checked by bibtool:

- Do not escape special characters like umlauts or accented characters. Instead, use the unicode character directly.
- You do not need to use braces to preserve capitalization as `DocumenterCitations.jl` keeps entries as is (in contrast to `bibtex`). In some cases, braces can even be harmful, i.e., show up in the output.
- If a DOI is available for your reference, please add it as a `doi` field to the BibTeX entry. In this case, please refrain from adding an additional `url` field.
- If your reference has no DOI or the paper is not open-access, but is available as an arXiv preprint, you can add the arXiv link as a `eprint` field (even additionally to a `doi` field). For other preprint servers (e.g. HAL), please refer to the [DocumenterCitations.jl docs](https://juliadocs.org/DocumenterCitations.jl/stable/syntax/#Preprint-support).
- Documents available only as an arXiv preprint should be added as `@Misc` entries with the arXiv-ID in the `eprint` field, e.g., `archiveprefix = {arXiv}` and `eprint = {2008.12651}`.
