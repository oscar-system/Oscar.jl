# Introduction for new developers

This document is meant to get new developers started. In general you have to do
the following six steps for submitting changes to the Oscar source:

1. Fork the main Oscar.jl repository. For this go to the [Oscar.jl github
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
4. Edit your source and try out your changes locally. To use your local copy of
   the sources, start Julia and
   ```
   ]dev /path/to/local/clone/of/your/fork/of/Oscar.jl
   ```
   If this succeeds, you can `using Oscar` in Julia and it will use your local
   copy.
5. Once you are done editing, push your branch and open a pull request. It is
   recommended that you open a draft pull request to the main Oscar repository
   as soon as you start working. That way Oscar developers are aware of work
   being done and can give feedback early in the process.
6. Once you have finished your work, mark your pull request as ready. It will
   then be reviewed and, probably after feedback and requests for changes,
   merged.


## The edit process
### Editing the source
The sources can be found in the `src` folder. Please pay attention to the
folder structure and choose sensibly where to place your code (when fixing a
bug this is probably a minor question).

### Adding tests

### Adding documentation
There are two places where documentation can added:
1. In the docstrings above the functions in the `src` folder;
2. In the documentation files in the `docs/src` folder.
In general, 1 is preferred to 2, i.e. any explanation of the functions and
objects should go there and the files in `docs/src` should remain relatively
sparse. Please also pay attention to the documentation section of [our style
guide](@ref styleguide).


## Further hints

### Ask Oscar related questions in the Oscar slack

### Use `]up`
Working with the development version also entails that the packages Oscar
depends on need to be up to date. Julia can update these packages if you type
`]up` in the Julia prompt. Many error messages after updating the source can be
resolved by simply updating.

### Style guide
Please have a look at [our style guide](@ref styleguide) to get an overview
over naming conventions, code formatting, etc.

### Building the documentation
To build and test the documentation, please have a look at [this page](@ref
building_docs).

### Rebasing
One way to stay up to date with the current master is rebasing. In order to do
this, add the main Oscar.jl repository as a remote, fetch, and then rebase.
```
git remote add oscar-system git@github.com:oscar-system/Oscar.jl
git fetch oscar-system
git rebase oscar-system/master
```
Adding the remote only has to be executed once.
