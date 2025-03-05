# Release management

We roughly follow the julia release process. The following outlines the approach taken for OSCAR 1.0.0 but is subject to change depending in how this works out.

Before preparing the release branch two developers are chosen to take care of managing the branches and taking final decisions on which changes can be merged during the release process.

## Branches

Assuming the current master has version number `X.Y.0-DEV`, in preparation for the next release
a branch `release-X.Y` is created directly on GitHub.
Once that is done the master should be set to `X.Y+1.0-DEV`.

!!! note 
    Please use the script `etc/update_version.sh` for version changes.

## Backports

The backports should be mostly bug-fixes and small improvements but no major new features or refactorings, to minimize the risk of regressions.
Every change that is considered for backporting has to be merged to master first.
Once in the stage of release candidates only small bug-fixes and documentation changes should be added.

Changes for the release, after the branch was created, are mostly done by cherry-picking the commits from master that correspond to the already _merged_ PRs. Preferably with `git cherry-pick -x` to keep a reference to the original commit.
The cherry-picked commits are collected on a `backports-release-X.Y` branch and a corresponding PR, like for [release 1.0.0-rc1](https://github.com/oscar-system/Oscar.jl/pull/3378). Please take care to set the merge target for the PR to the `release-X.Y` branch.

Once this set of changes is complete the PR is merged to the release branch, using a _proper merge commit_, i.e. please do not squash or rebase this.
The list of backported PRs may be added to the commit message for the merge.

### Conflicts

If changes from another PR cannot be cleanly applied to the backports branch it may be necessary to create a separate PR targeting the _backports branch_ to manually backport changes from that PR and resolve any conflicts. Make sure that the branch for this PR branches off from the release branch and not from master (to avoid picking up unrelated commits). Then cherry-pick the relevant commits and fix any conflicts.

## Pre-releases

For the next pre-release create a new PR to the release branch adjusting the version, this may be merged with a rebase to avoid some extra commits. That commit can then be tagged on GitHub as a new pre-release.

!!! note
    Pre-releases cannot be registered in the julia registry.

For OSCAR 1.0.0 we plan to have one or two release candidates but unlike julia and due to the tight schedule no `-alpha` or `-beta` versions.

## Releases

Just before making the release, `CHANGELOG.md` should be updated.
See [Updating `CHANGELOG.md`](@ref) for details.

Once the version is at X.Y.0, the version can be registered in the julia registry with `@JuliaRegistrator register()`. In this case TagBot will create the corresponding release, but preferably recheck and clean up the list of changes since unfortunately TagBot currently does not take the branch into account and will show all PRs that have been merged to master.

Consider deploying a released version at [zenodo](https://zenodo.org/).

## Bugfix releases

Further bugfix releases X.Y.Z can be created in a similar manner, skipping the pre-releases.
