# Serialization

!!! warning
    Never load data from an untrusted source. Loading data is inherently unsafe
    and at this point allows arbitrary code execution on your machine. Just as
    you should never run a program from someone you do not trust, you should
    also not load their data.

!!! warning
    Serialization development has just started and the concrete design may
    still change drastically. In particular the mechanism for upgrading old
    data to newer versions is not in place yet, so at this point we do not yet
    guarantee that old data can be read.

This document summarises the serialization efforts of OSCAR, how it is supposed
to work, how it works and the overall goal.
[Serialization](https://en.wikipedia.org/wiki/Serialization) broadly speaking
is the process of writing data to and reading data from files. There are many
reasons for needing this feature in OSCAR, but the main reason is communication
on mathematics by mathematicians.

## How it works
The mechanism for saving and loading is very simple. It is implemented via two
methods `save` and `load`, and works in the following manner:
```
julia> save("/tmp/fourtitwo.json", 42);

julia> load("/tmp/fourtitwo.json")
42

```
As hinted by the filename, OSCAR writes a file in JSON format. The file looks
as follow:
```
{
    "_ns": {
        "Oscar": [
            "https://github.com/oscar-system/Oscar.jl",
            {
                "major": 0,
                "minor": 8,
                "patch": 3,
                "prerelease": [
                    "DEV"
                ],
                "build": []
            }
        ]
    },
    "id": "-1",
    "type": "Base.Int",
    "data": "42"
}
```
It contains the version of OSCAR it was written by, its type, and the actual
content, in this case as a string.

### The `id`
If you look at the file `src/Serialization/main.jl`, you will see that all
`save` methods hand down a `SerializerState`, and all `load` methods have a
`DeserializerState`. These two objects are very simple, they just contain
dictionaries connecting objects and their `id`. We use this to avoid saving or
loading larger objects twice (or multiple times). Consider the following
example code snippet:
```
c = cube(3);
LP0 = LinearProgram(c, [2,2,-3]);
LP1 = LinearProgram(c, [2,2,4]);
v = [LP0, LP1];
save("vector_of_lp.json", v)
```
This creates two linear programs on the cube, stores them in a vector and then
writes this vector to a file. It would be wasteful to store the cube twice for
each linear program, instead it is only stored once and the second linear
program just gets the `id` of the cube in its serialized form. Please take some
time to look at the file written in this concrete example.

### The version number
We will use the version number for checking compatibility of the data with the
current OSCAR version before attempting to load it. If the data version is
lower than the OSCAR version we will provide appropriate upgrade scripts such
that the data can be loaded. We will not provide scripts for attempting to
downgrade data, but we will throw a warning or even error in this case. We may
provide an option for attempting to load anyway in such a scenario.

### Implementation
All files for serialization can be found in the folder `src/Serialization`. The
naming conventions of the files there follows the overall structure of OSCAR,
i.e. the file `src/Serialization/PolyhedralGeometry.jl` contains functions for
serializing objects of the polyhedral geometry section.

The file `main.jl` contains the core of the serialization process, namely:
- reading and writing files;
- the `SerializerState` and `DeserializerState` objects;
- writing and reading versions; and
- generic functions for attempting to serialize objects that do not have their
  own dedicated serialization methods.

If you want to write a serialization routine for an object, the way to go is to
implement the following two functions, here in the example for `fmpz`:
```
function load_internal(s::DeserializerState, ::Type{fmpz}, str::String)
    return fmpz(str)
end

function save_internal(s::SerializerState, z::fmpz)
    return string(z)
end
```
Then the main serialization methods will dispatch to `load_internal` and
`save_internal` for `fmpz` instead of attempting the generic serialization.

Often the generic serialization will fail and it is necessary to provide a
`save_internal` and `load_internal` function. In that case, please have a look at
the existing functions to get an idea of how these work, and maybe use
something of this as a blueprint.

## Challenges
This section documents the various challenges we (will) encounter while
implementing this feature.
- OSCAR is based on several subsystems, some of which already have their own
  serialization. We want this to be compatible, if possible in both directions.
- Many mathematical objects need context to be understood. A polynomial needs
  the ring it lives in, a group element needs the surrounding group, a divisor
  needs the underlying variety, etc. We will need a way to store this context
  along the objects.
- Context should not be stored twice: A matrix of polynomials should only store
  the surrounding ring once.
- Support other data formats: It has been proposed to not only support JSON,
  but binary formats needed for HPC communication as well. It is unclear
  whether this needs a separate implementation.
- Versioning and upgrading: Work on OSCAR will change what its objects look
  like. Nevertheless, we still want to be able load data written by older
  versions of OSCAR. For this we intend to develop an upgrade mechanism.

Another important point is the wider mathematical context of the data and code.
For data associated to a publication, this context is provided by the paper.



## Goals

The general goal is to make mathematical data
[FAIR](https://en.wikipedia.org/wiki/FAIR_data), a goal for which we cooperate
with the [MaRDI](https://www.mardi4nfdi.de/about/mission) project.

The ramifications of making mathematical data FAIR are manifold. 
- It becomes easier to exchange data and code with fellow mathematicians,
  enhancing communication and boosting research.
- Computer experiments and new implementations require a lot of work and hence
  deserve to be recognized in form of a publication. Standardizing data plays
  an important role for this process.
- Future generations of mathematicians will be able to reuse both data and code
  if we establish a FAIR culture.
