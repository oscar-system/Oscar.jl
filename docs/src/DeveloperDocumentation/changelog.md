```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Updating `CHANGELOG.md`

Before every release, `CHANGELOG.md` should be updated. Ideally most of this is taken care of
by automation, provided PR labels were sensibly set. Below we explain which
labels we use; how people who create or merge PRs should set the labels; and how to
perform the actual update.

## How to set labels when submitting or reviewing a PR

Everyone submitting a PR should ideally follow this checklist if they also have the access
rights for assigning labels. If they don't, they should ideally leave instructions that
indicate which labels to set, and a PR reviewer should then perform the following steps.

1. If your PR is *not* supposed to appear in the changelog, assign the `release notes: not needed` label and stop here.
2. If possible, edit the title of your PR to be usable as an entry in the changelog, and assign the `release notes: use title` label.
   - Hint: Try to start the title with a verb in present tense such as "Add", "Enhance", "Remove".
   - Hint: If function or type names appear in the title, surround them by backticks.
3. If you need a longer description than fits into the title, leave a comment with your suggestion
   for this description, and assign the `release notes: to be added` label.
4. In either case, also assign at least one `topic: XYZ` label (see below for an overview),
   and possibly some of the other labels described in the following section as you deem fit.

Optionally, if you think your change is exceptionally noteworthy, also assign the
`release notes: highlight` label.

### Example

For example, if the change is a PR with title `Bump GAP.jl to 0.13`, and labels `package: GAP` and
`release notes: use title`, it would show up in the changelog file as

    ### Changes related to the package GAP

    - [#4421](https://github.com/oscar-system/Oscar.jl/pull/4421) Bump GAP.jl to 0.13


## How we use pull request labels to organize the changelog

The script `dev/releases/release_notes.py` will generate updates to `CHANGELOG.md`
automatically, based on PRs merged after the previous release, and what those PRs are
labelled. We have the following labels, along with how they are meant to be applied.

### Primary Labels: Release Notes Behaviour

| Label | Meaning |
|-------|---------|
| `release notes: added`        | The release notes for this PR were manually added to the changelog, and should be ignored by the script |
| `release notes: not needed`   | This PR does not warrant an entry in the release notes. Internal only changes, like reorganization of private functions, changes to the test pipeline, etc can be tagged with this |
| `release notes: use title`    | The release notes for this PR should be based on the title of this PR. The script will turn \$TITLE from the PR to `[#xyz] $TITLE` |
| `release notes: to be added`  | These PRs will be marked by the script as a prominent TODO item. Check these PRs manually, and after updating them appropriately, relabel these items to either `release notes: added` or `release notes: use title` |
| none of the above             | These PRs will be added to a separate prominent TODO category. Check these PRs manually, and after updating them appropriately, relabel these items to one of `release notes: added`, `release notes: use title`, or `release notes: not needed` |

### Secondary Labels: Topic

In addition to the release notes action labels, you can tag your PR with these following
labels, and the release notes script will organize them appropriately:

| Label                         | Changelog Category |
|-------------------------------|--------------------|
| `release notes: highlight`    | Highlights |
| `renaming`                    | Items being renamed |
| `topic: algebraic geometry`   | Changes related to Algebraic Geometry |
| `topic: combinatorics`        | Changes related to Combinatorics |
| `topic: commutative algebra`  | Changes related to Commutative Algebra |
| `topic: FTheoryTools`         | Changes related to F-Theory Tools |
| `topic: groups`               | Changes related to Groups |
| `topic: lie theory`           | Changes related to Lie Theory |
| `topic: number theory`        | Changes related to Number Theory |
| `topic: polyhedral geometry`  | Changes related to Polyhedral Geometry |
| `topic: toric geometry `      | Changes related to Toric Geometry |
| `topic: tropical geometry`    | Changes related to Tropical Geometry |
| `serialization`               | Changes related to serializing data in the MRDI file format ? |
| `enhancement`                 | New features or extended functionality |
| `experimental`                | Only changes experimental parts of OSCAR |
| `optimization`                | Performance improvements or improved testing |
| `bug: crash`                  | Fixed bugs that could lead to crashes |
| `bug`                         | Other fixed bugs |
| `documentation`               | Improvements or additions to documentation |
| `package: AbstractAlgebra`    | Changes related to the package AbstractAlgebra |
| `package: AlgebraicSolving`   | Changes related to the package AlgebraicSolving |
| `package: GAP`                | Changes related to the package GAP |
| `package: Hecke`              | Changes related to the package Hecke |
| `package: Nemo`               | Changes related to the package Nemo |
| `package: Polymake`           | Changes related to the package Polymake |
| `package: Singular`           | Changes related to the package Singular |

## Suggestions for formulations

In general the description of each change should start with a verb in present
tense. Here are some more concrete suggestions.

| Change                        | Example |
|-------------------------------|--------------------|
| move from experimental to src/ | Graduate bla from experimental to officially supported |
| feature added                 | Add `bla` for `blub `/ Support `bla` for `blub` / Implement `bla`
| renaming things               | Rename `bla` to `blub`
| bug fix                       | Fix `bla` in `blub`
| improvements                  | Improve (performance) of `blub`
| experimental feature          | Experimental: add support for `bla`

## Updating the changelog

There are two ways to update the changelog: by invoking the script directly, or by triggering
a GitHub workflow.

This is typically done by the release managers, i.e., you don't have to do this when
submitting a pull request. In particular, don't manually add anything to `CHANGELOG.md`!

### Changelog Script

The script is located at `dev/releases/release_notes.py`. It must be run from
inside the `dev/releases` directory. Running it updates `CHANGELOG.md`, for review.

### GitHub Workflow

There is a GitHub workflow which runs the `release_notes.py` script and turns the result
into a pull request. Anyone with write permission to the OSCAR repository can use it by
visiting <https://github.com/oscar-system/Oscar.jl/actions/workflows/changelog.yml>, then
selecting the "Run workflow" button. This should open a little drop down, offering to
select a branch (usually this should be left at the default `master`). Click the green
"Run worklow" button to start the process.

After running a bit, this will either open a new PR against `master`, or update an
existing PR created previously by it. This suggests the following usage pattern: run the
workflow to create an initial PR, and inspect it. If the result is not satisfactory (e.g.
because it contains uncategorized entries, or entries with typos, unclear language,
formatting issues, etc.), one can simply adjust the titles of those PRs, and then re-run
the workflow. Repeat as needed.
