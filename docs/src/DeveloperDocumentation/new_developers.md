# Introduction for new developers
This document is meant to get new developers started. It will not go into depth
of programming in Julia or working with git, as there are far better resources
on these things online.

!!! note "Pay attention to your GitHub notifications!"
    Once you open a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) on GitHub you will receive feedback, comments,
    and questions on GitHub. So please pay attention to your GitHub
    notifications.

## Important notes
1. If you encounter error messages after rebasing to the current master, chances
   are that some dependencies need upgrading. Please first try whether
   executing `]up` gets rid of your errors.
2. Please have a look at the [Developer Style Guide](@ref) and the [Design
   Decisions](@ref). Adhering to the style guide makes reviewing code easier
   for us, and hence your new feature can be merged faster.
3. Let us know what you are working on early:
   - You can open a draft [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) on GitHub right at the beginning of your
     work. We are more than happy to look at incomplete prototypes to get an
     idea of what you are working on. This allows us to assess what kind of
     problems you might encounter and whether we can mitigate these by making
     changes to OSCAR.
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
2. [Clone your forked repository to your local
   machine](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).
   If you have set up ssh access you can do this in the following way:
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
5. Once you are done editing, push your branch and open a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request). It is
   recommended that you open a draft [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to the main Oscar repository
   as soon as you start working. That way Oscar developers are aware of work
   being done and can give feedback early in the process.
6. Once you have finished your work, mark your [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) as ready. It will
   then be reviewed and, probably after feedback and requests for changes,
   merged.

## Alternative: `]dev Oscar`

Alternatively you can call
```
]dev Oscar
```
in Julia. This will create a directory `~/.julia/dev/Oscar`. This directory is
a git clone of the central Oscar repository. You can develop your code here,
however you will still have to fork Oscar, as you have no rights to push to the
central repository. You can then add your fork as another remote, have a look
at the section on rebasing below for hints.


## The edit process
### Editing the source
The sources can be found in the `src` folder. Please pay attention to the
folder structure and choose sensibly where to place your code (when fixing a
bug this is probably a minor question).

### Adding tests
There are two ways to add tests:
  - There are combined tests and examples in the docstrings, namely the
    `jldoctest` blocks. For these have a closer look at [Documenting OSCAR
    code](@ref).
  - Larger tests and tests that aren't useful examples are in the `test`
    folder. The main file there is `test/runtests.jl` which then includes other
    testfiles. 

```@docs
Oscar.test_module
```

### Adding documentation
For more information on docstrings, please read our page on [Documenting OSCAR
code](@ref).  There are two places where documentation can be added:
1. In the docstrings above the functions in the `src` folder;
2. In the documentation files in the `docs/src` folder. The overall structure
   is fixed in the file `docs/doc.main`. If you create a new file in
   `docs/src`, you will have to add an entry in `docs/doc.main`.

In general, 1 is preferred to 2, i.e. any explanation of the functions and
objects should go there and the files in `docs/src` should remain relatively
sparse. Please also pay attention to the documentation section of the
[Developer Style Guide](@ref).


## Further hints

### Give [gh](https://github.com/cli/cli) a try
Especially if you will be developing a lot, this can speed up your workflow
tremendously.

### Use the [`Revise`](https://github.com/timholy/Revise.jl) package
Using `Revise` you can avoid having to restart Julia and reloading OSCAR when
editing the code. As a quick summary, first install `Revise` with:
```
using Pkg;
Pkg.add("Revise");
```
From then on do
```
using Revise,Oscar
```
whenever you are using OSCAR in Julia.

### Ask Oscar Related Questions in the [Oscar Slack](https://join.slack.com/t/oscar-system/shared_invite/zt-thtcv97k-2678bKQ~RpR~5gZszDcISw).

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
