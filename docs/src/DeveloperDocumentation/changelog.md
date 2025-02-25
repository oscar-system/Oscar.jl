# Updating `CHANGELOG.md`

Before every release, `CHANGELOG.md` should be updated. Most of this is ideally
taken care by automation, provided PR labels were sensibly set. The following
labels explain which labels we use; how people who create or merge PRs should
set the labels; and how to perform the actual update.

## How we use pull request labels to organize the changelog

The script [release_notes.py](../../../dev/releases/release_notes.py) will generate updates to 
`CHANGELOG.md` automatically, based on PRs merged after the previous release, and what those PRs
are labelled. We have the following labels, along with how they are meant to be applied

### Primary Labels: Release Notes Behaviour

| Label | Meaning |
|-------|---------|
| release notes: added | The release notes for this PR were manually added to the changelog, and should be ignored by the script |
| release notes: not needed | This PR does not warrant an entry in the release notes. Internal only changes, like reorganization of private functions, changes to the test pipeline, etc can be tagged with this |
| release notes: use title | The release notes for this PR should be based on the title of this PR. The script will turn $TITLE from the PR to `[#xyz] $TITLE` |
| release notes: to be added | These PRs will be marked by the script as a prominent TODO item. Check these PRs manually, and after updating them appropriately, relable these items to either `release notes: added` or `release notes: use title` |
| \<no label\> | These PRs will be added to a separate prominent TODO category. Check these PRs manually, and after updateing them appropriately, relable these items to one of `release notes: added`, `release notes: use title`, or `release notes: not needed` |

### Secondary Labels: Topic

In addition to the release notes action labels, you can tag your PR with these following labels, and the release notes script will organize them appropriately :

| Label | Changelog Category |
|-------|--------------------|
| enhancement | New features or extended functionality |
| experimental | Only changes experimental parts of OSCAR |
| optimization | Performance improvements or improved testing |
| package: AbstractAlgebra | Changes related to the package AbstractAlgebra |
| package: AlgebraicSolving | Changes related to the package AlgebraicSolving |
| package: GAP | Changes related to the package GAP |
| package: Hecke | Changes related to the package Hecke |
| package: Nemo | Changes related to the package Nemo |
| package: Polymake | Changes related to the package Polymake |
| package: Singular | Changes related to the package Singular |
| renaming | Items being renamed |
| serialization | Changes related to serializing data in the MRDI file format ? |
| topic: algebraic geometry | Changes related to Algebraic Geometry |
| topic: combinatorics | Changes related to Combinatorics |
| topic: FTheoryTools | Changes related to F-Theory Tools |
| topic: groups | Changes related to Groups |
| topic: LieTheory | Changes related to Lie Theory |
| topic: number theory | Changes related to Number Theory |
| topic: polyhedral geometry | Changes related to Polyhedral Geometry |
| topic: rings | Changes related to Rings |
| topic: schemes | Changes related to Schemes |
| topic: toric schemes | Changes related to Toric Schemes |
| topic: toric varieties | Changes related to Toric Varieties |
| topic: tropical geometry | Changes related to Tropical Geometry |
| bug: crash | Fixed bugs that could lead to crashes |
| bug | Other fixed bugs |
| documentation | Improvements or additions to documentation |

## Example

For example, if you have a PR titled `Bump GAP.jl to 0.13`, and labelled it `package: GAP`, and
`release notes: use title`, it would show up in the changelog file as


    ### Changes related to the package GAP
    
    - [#4421](https://github.com/oscar-system/Oscar.jl/pull/4421) Bump GAP.jl to 0.13



## Suggestions for categorizing pull requests

todo ? Or is the above enough ?

## Updating the changelog

There are two ways update the changelog: invoking the script directly, or calling the github action.

### Changelog Script

The script is located at
[../../../dev/releases/release_notes.py](../../../dev/releases/release_notes.py). Running it will
update `CHANGELOG.md`, for review.

### Github Workflow

There is a github CI workflow at
https://github.com/oscar-system/Oscar.jl/actions/workflows/changelog.yml . This will run the script,
and automatically make a PR against `master`, updating the changelog. If there is an already open
PR by the workflow, running it again will update this PR.
