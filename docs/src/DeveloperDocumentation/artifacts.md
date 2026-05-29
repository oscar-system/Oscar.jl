```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Artifacts

This page explains what artifacts are, when to use them, and how to work with them in `Oscar.jl`.



## What Artifacts Are and When to Use Them

Artifacts are content-addressed bundles managed by [Julia's artifact system](https://pkgdocs.julialang.org/v1/artifacts/) and declared in `Oscar.jl/Artifacts.toml`. Artifacts are automatically downloaded, or updated when `Oscar.jl` is installed or updated, unless the required artifact is already present locally.

Julia also supports [lazy artifacts](https://docs.julialang.org/en/v1/stdlib/LazyArtifacts/), which are installed only on demand. In the following, we shall ignore lazy artifacts.

Artifacts allow to:

- keep the `Oscar.jl` repository small,
- distribute data reproducibly,
- avoid duplication of large files across versions.

Artifacts should be used for data that does not belong directly in the `Oscar.jl` repository, in particular:

- generated data,
- large datasets,
- data not intended to be edited manually.

As a rule of thumb, data stored directly in the `Oscar.jl` repository should not exceed the size of typical source files (around 100 KB). Small examples, test inputs, and simple hand-written data should remain in the `Oscar.jl` repository.



## Creating, Hosting and Using Artifacts

This section describes how to create, host, register, and use artifacts in the `Oscar.jl` repository. The workflow typically proceeds as follows:

- [ ] [Serialize data](@ref data_preparation)
- [ ] [Upload the tarball to a stable hosting location](@ref artifact_hosting)
- [ ] [Add an entry to `Oscar.jl/Artifacts.toml` and open a pull request](@ref artifact_registration)
- [ ] [Use the artifact once the pull request has been merged](@ref artifact_usage)

We illustrate this workflow with an [end-to-end example](@ref end-to-end-artifact-example).


### [End-to-End Example](@id end-to-end-artifact-example)

The following example illustrates the complete workflow for creating, registering, and using an artifact in `Oscar.jl`.

#### Create data

```julia
using Oscar

obj = ZZ(2)^10
```

#### Save data

Let us employ the  [`.mrdi` file format to serialize our data](https://docs.oscar-system.org/dev/DeveloperDocumentation/serialization/). (You may of course use the file format you prefer.)

```julia
save("example.mrdi", obj)
```

#### Pack into a tarball

Create a compressed tarball, for example `example_data_v1.tar.gz`, containing the file `example.mrdi`.

!!! note
    Do not introduce an intermediate directory inside the tarball. If multiple files are included, place them directly in the archive.

#### Register the artifact

Register the artifact by adding an entry to `Oscar.jl/Artifacts.toml`; see also [Julia's artifact documentation](https://pkgdocs.julialang.org/v1/artifacts/). The package [`ArtifactUtils.jl`](https://github.com/JuliaPackaging/ArtifactUtils.jl) can help automate parts of this workflow. The following text demonstrates the manual workflow.

First, compute the `sha256` and the `git-tree-sha1`:

```julia
using Tar, Inflate, SHA

filename = "/absolute/path/to/example_data_v1.tar.gz"

println("sha256: ", bytes2hex(open(sha256, filename)))
println("git-tree-sha1: ", Tar.tree_hash(IOBuffer(inflate_gzip(filename))))
```

Then add this information, together with the host location, to `Oscar.jl/Artifacts.toml`. This very example, we host at [https://martinbies.github.io/Materials/Data/example_data_v1.tar.gz](https://martinbies.github.io/Materials/Data/example_data_v1.tar.gz). Then, the corresponding entry to `Oscar.jl/Artifacts.toml` takes the following form:

```toml
[MyExample]
git-tree-sha1 = "ff2f21e623a130f47116d847ae54fd55232b42c1"
lazy = true

    [[MyExample.download]]
    sha256 = "3bd5a20e84b2e579ecde7dc5d7d4606444daf08407eecc9d8e59be1e468ca5a1"
    url = "https://martinbies.github.io/Materials/Data/example_data_v1.tar.gz"
```

You may of course replace `MyExample` in the above entry with any other string that you find descriptive.

#### Use the artifact

The artifact string macro is exported by `LazyArtifacts`. In a standalone Julia session, load it explicitly.

```julia
using LazyArtifacts

obj_path = artifact"MyExample/example.mrdi"

obj = load(obj_path)
```

In `Oscar.jl` source files, the artifact string macro is already available and `using LazyArtifacts` is typically not required.


### [Preparing Data (Serialization)](@id data_preparation)

Creating artifacts requires that the corresponding data be serialized locally first. Details are provided in the [Serialization page](https://docs.oscar-system.org/dev/DeveloperDocumentation/serialization/).

We recommend the use of the `.mrdi` file format for serialization. However, this is not a strict requirement and you may use any file format that you see fit.


### [Hosting Artifacts](@id artifact_hosting)

Artifacts should be hosted at stable and persistent locations.

Preferred options include:

- archival services such as [Zenodo](https://zenodo.org/communities/oscar/records?q=&l=list&p=1&s=10&sort=newest), in particular for long-term or publication-related data,
- other stable hosting solutions agreed upon by the maintainers.

For historical reasons, some artifacts are currently hosted via GitHub release assets, for example at
[Oscar.jl/archive-tag-1](https://github.com/oscar-system/Oscar.jl/releases/tag/archive-tag-1). This approach should be used with care, as GitHub release assets are not intended to function as a long-term artifact registry.

!!! warning
    Ensure that existing files are never removed or renamed, as they may be required by older `Oscar.jl` releases.

When debugging, contributors are encouraged to use temporary staging areas before publishing long-term artifact versions. In particular, publication-related Zenodo entries should typically not be cluttered with intermediate or broken artifact versions created during development.

Data intended for querying may also be suitable for [OscarDB](https://docs.oscar-system.org/dev/Experimental/OscarDB/introduction/).


### [Registering](@id artifact_registration)

Recall that creating an artifact typically involves the following steps:

- collecting the relevant data files,
- packing these files into a compressed tarball (`.tar.gz`),
- uploading the tarball to a stable hosting location,
- adding a corresponding entry to `Oscar.jl/Artifacts.toml`,
- opening a pull request with the change to `Oscar.jl/Artifacts.toml`.

Registering refers to the final two steps. Once the change to `Oscar.jl/Artifacts.toml` is merged into the `Oscar.jl` repository, the artifact becomes available.

The [end-to-end example](@ref end-to-end-artifact-example) explicitly demonstrates the required changes to `Oscar.jl/Artifacts.toml`. Additional information is available in [Julia's artifact documentation](https://pkgdocs.julialang.org/v1/artifacts/).


### [Using Artifacts in Oscar](@id artifact_usage)

Artifacts can be accessed via Julia's artifact system through the `artifact"..."` string macro. For example:

```julia
using LazyArtifacts

model_data_path = artifact"FTM-1511-03209/1511-03209.mrdi"

model = load(model_data_path)
```

The string passed to `artifact"..."` is determined by the corresponding entry in `Oscar.jl/Artifacts.toml`.

In `Oscar.jl` source files, the artifact string macro is already available and `using LazyArtifacts` is typically not required.



## Updating Artifacts

### General Rules

Artifacts are immutable. Updating an artifact therefore requires:

- creating a new tarball with the updated data,
- uploading the tarball to a stable hosting location,
- recomputing the `sha256` and `git-tree-sha1`,
- updating the corresponding entry in `Oscar.jl/Artifacts.toml`,
- opening a pull request with the changes to `Oscar.jl/Artifacts.toml`.

Once the pull request is merged, the updated artifact becomes available.

Updating an artifact on a hosting platform alone, for example by uploading a new version to Zenodo, is not sufficient. Any change to the artifact contents changes its hashes and therefore requires a corresponding update of `Oscar.jl/Artifacts.toml`.

!!! warning
    Any files referenced by the `Oscar.jl` master branch must not be modified, renamed, or deleted, as they may be required by earlier `Oscar.jl` releases.

When repeatedly debugging or refining artifacts, contributors are encouraged to use temporary staging areas before publishing long-term versions. In particular, publication-related Zenodo entries should typically not be cluttered with intermediate or broken artifact versions created during development.

### Serialization Upgrades

Note that the chosen file format for serialization may be subject to development. In particular, the `.mrdi` file format, which we recommend for serialization, is under active development. Consequently, its standard evolves over time. Older files remain compatible with newer OSCAR versions; however, loading older artifacts may require upgrade steps during deserialization. For large artifacts, these upgrades may become time consuming. It is therefore recommended to use the most recent serialization standard when creating artifacts and to periodically upgrade older artifacts if appropriate.

Additional details are provided in the [Serialization documentation](https://docs.oscar-system.org/stable/DeveloperDocumentation/serialization/#Upgrades).



## Artifacts in Upstream Dependencies of `Oscar.jl`

### GAP.jl

`GAP.jl` uses artifacts to install [GAP](https://www.gap-system.org/) and GAP packages. Detailed maintainer information is provided in [`GAP.jl/README.maintainer.md`](https://github.com/oscar-system/GAP.jl/blob/master/README.maintainer.md).
