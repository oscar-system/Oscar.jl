```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# [Serialization](@id dev_serialization)
This document summarizes the serialization efforts of OSCAR, how it works, and what our long-term vision is.
[Serialization](https://en.wikipedia.org/wiki/Serialization) broadly speaking
is the process of reading and writing data. There are many reasons for this
feature in OSCAR, but the main reason is communication on mathematics by
mathematicians.

We implement our serialization in accordance with the [MaRDI](https://www.mardi4nfdi.de/about/mission) file format specification developed by Della Vecchia, Joswig and Lorenz [D-VJL24*1](@cite).
In particular, we use a JSON extension to serialize data.


## How it works
The mechanism for saving and loading is very simple. It is implemented via two
methods `save` and `load`, and works in the following manner:
```julia-repl
julia> save("/tmp/fourtitwo.mrdi", 42);

julia> load("/tmp/fourtitwo.mrdi")
42

```
The filename hints to the MaRDI file format [D-VJL24*1](@cite), which employs JSON.  The file looks as follows:
```json
{
  "_ns": {
    "Oscar": [
      "https://github.com/oscar-system/Oscar.jl",
      "0.14.0-DEV-8fe2abbe39890a7d3324adcba7f91812119c586a"
    ]
  },
  "_type": "Base.Int",
  "data": "42"
}
```
It contains the precise version of OSCAR used for this serialization.
The content is "42", it represents a `Base.Int`, according to the `_type` field.


## Implementation
To list and describe all implementations and encodings of all types in OSCAR
is not a possible feat due to the arbitrarily deep and nested type structures
available in OSCAR, we point any developer looking to understand the encodings
of certain types to the OSCAR source code. All files for serialization can be
found in the folder `src/Serialization`. The convention of the files there
follows the overall structure of OSCAR, i.e. the file
`src/Serialization/PolyhedralGeometry.jl` contains functions for
serializing objects of the polyhedral geometry section.

We include a basic example of the encoding for a `QQPolyRingElem` so one can get a taste
before delving into the source code. Here we store the polynomial `x^3 + 2x + 1//2`.
The encoding for polynomials is to store a list of tuples, where each entry in the
list represents a term of the polynomial and where the first entry of the tuple is
the exponent and the second entry is the coefficient. Here we serialize a univariate
polynomial so the first entries are always integers, in general this may be an array
of integers. The coefficients here are elements of `QQ` however in general
the coefficients themselves may be described as polynomials of the generators
of some field extension, i.e. the second entry may again be a list of tuples and so on.
The nested structure of the coefficient will depend on the description of the field
extension.


```json
{
  "_ns": {
    "Oscar": [
      "https://github.com/oscar-system/Oscar.jl",
      "1.1.0-DEV-6f7e717c759f5fc281b64f665c28f58578013c21"
    ]
  },
  "_refs": {
    "e6c5972c-4052-4408-a408-0f4f11f21e49": {
      "_type": "PolyRing",
      "data": {
        "base_ring": {
          "_type": "QQField"
        },
        "symbols": [
          "x"
        ]
      }
    }
  },
  "_type": {
    "name": "PolyRingElem",
    "params": "e6c5972c-4052-4408-a408-0f4f11f21e49"
  },
  "data": [  [      "0",      "1//2"    ],
	     [      "1",      "2"       ],
	     [      "3",      "1"       ]  ]
}

```

When trying to understand an encoding of a particular type in OSCAR it is always
best to make a minimal example and store it. Then use a pretty printer to format
the JSON, and have it close by while going through the source code.

### Description of the saving and loading mechanisms

We require that any types serialized through OSCAR are registered using
[`@register_serialization_type`](@ref).
This is to ensure user safety during the load process by avoiding code
evaluation.

```@docs
@register_serialization_type
```

There are three pairs of saving and loading functions that are used
during serialization:
1. `save_typed_object`, `load_typed_object`
2. `save_type_params`, `load_type_params`
3. `save_object`, `load_object`


#### `save_type_object` / `load_type_object`

These functions should not be touched, they are high level
functions and are used to (de)serialize the object with its
type information as well as its data. The data and type nodes are
set in `save_typed_object` resulting in a "data branch" and "type branch".


#### `save_type_params` / `load_type_params`

The serialization mechanism stores data in the format of a tree, with the
exception that some nodes may point to a shared reference. The "data branch"
is anything that is a child node of a data node, whereas the "type branch" is
any information that is stored in a node that is a child of a type node.

These functions should also not be touched, however they expect an implementation
of `type_params` whenever saving a type `T`. By default `type_params` will return
`nothing`. The `type_params` function does a shallow pass through an `obj` of type
`T` gathering the necessary parameters for serializing `obj`.
In most cases these parameters are the parameters of the `obj` that uses references.
For example if `obj` is of type `RingElem` than it is expected that `type_params`
should contain at least 

#### `save_object` / `load_object`

These functions should be the first functions to be overloaded when
implementing the serialization of a new type.
The functions `save_data_dict` and `save_data_array` are helpers functions
that structure the serialization.

The examples show how they can be used to save data using the structure
of an array or dict. Each nested call to `save_data_dict` or `save_data_array`
should be called with a key that can be passed as the second parameter.

##### Examples

###### Example 1
```julia
function save_object(s::SerializerState, obj::NewType)
  save_data_array(s) do
    save_object(s, obj.1)
    save_object(s, obj.2)
   save_data_dict(s) do
      save_object(s, obj.3, :key1)
      save_object(s, obj.4, :key2)
    end
  end
end
```

This will result in a data format that looks like this.
```json
[
  obj.1,
  obj.2,
  {
    "key1": obj.3,
    "key2": obj.4
  }
]
```

With the corresponding loading function similar to this.

```julia
function load_object(s::DeserializerState, ::Type{<:NewType})
  (obj1, obj2, obj3_4) = load_array_node(s) do (i, entry)
    if entry isa JSON3.Object
      obj3 = load_object(s, Obj3Type, :key1)
      obj4 = load_object(s, Obj3Type, :key2)
      return OtherType(obj3, obj4)
    else
      if p(entry) == c
        load_object(s, Obj1Type)
      else
        load_object(s, Obj2Type)
      end
    end
  end
  return NewType(obj1, obj2, obj3_4)
end
```

##### Example 2
```julia
function save_object(s::SerializerState, obj::NewType)
  save_data_dict(s) do
    save_object(s, obj.1, :key1)
    save_data_array(s, :key2) do
      save_object(s, obj.3)
      save_typed_object(s, obj.4) # This is ok
    end
  end
end
```
This will result in a data format that looks like this.

```json
{
  "key1": obj.1,
  "key2":[
    obj.3,
    {
      "type": "Type of obj.4",
      "data": obj.4
    }
  ]
}
```

The corresponding loading function would look something like this.
```julia
function load_object(s::DeserializerState, ::Type{<:NewType}, params::ParamsObj)
   obj1 = load_object(s, Obj1Type, params[1], :key1)

   (obj3, obj4) = load_array_node(s, :key2) do (i, entry)
     if i == 1
       load_object(s, Obj3Type, params[2])
     else
       load_typed_object(s)
     end
   end
   return NewType(obj1, OtherType(obj3, obj4))
 end
```

This is ok
```julia
function save_object(s::SerializerState, obj:NewType)
  save_object(s, obj.1)
end
```

While this will throw an error
```julia
function save_object(s::SerializerState, obj:NewType)
  save_object(s, obj.1, :key)
end
```

If you insist on having a key you should use a `save_data_dict`.
```julia
function save_object(s::SerializerState, obj:NewType)
  save_data_dict(s) do
    save_object(s, obj.1, :key)
  end
end

function load_object(s::SerializerState, ::Type{<:NewType})
  load_node(s, :key) do x
    info = do_something(x)

    if info
      load_object(s, OtherType)
    else
      load_object(s, AnotherType)
    end
  end
end
```

Note for now `save_typed_object` must be wrapped in either a `save_data_array` or
`save_data_dict`. Otherwise you will get a key override error.

### Serializers

The code for the different types of serializers and their states is found in the
`serializers.jl` file. Different serializers have different use cases, the
default serializer `JSONSerializer` is used for writing to a file.

When passing `serialize_refs = false` to the `JSONSerializer` it will not
store the refs of any types that are registered with the `uses_id` flag.
When using this flag it is left up to the user to guarantee that any refs
mentioned in the file are loaded prior to loading the file.
This is useful for cases where the user wants to store multiple objects that
refer to the same object, but does not want to store the refs in each file.
Instead, one can now store the refs in a separate file, and store the objects
themselves without the refs.

There is also the `IPCSerializer` which at the moment is
equal to the `JSONSerializer(serialize_refs = false)`. However, this serializer
may be changed in the future to support binary representations of some types
for faster inter-process communication (IPC).

### Upgrades

All upgrade scripts can be found in the `src/Serialization/Upgrades` folder.
The mechanics of upgrading are found in the `main.jl` file where the
[`Oscar.Serialization.upgrade`](@ref) function provides the core functionality. Upgrading
is triggered during [`load`](@ref) when the version of the file format
to be loaded is older than the current OSCAR version.

```@docs
Oscar.Serialization.upgrade
Oscar.Serialization.upgrade_data
```

#### Upgrade Scripts

All upgrade scripts should be contained in a file named after the version
they upgrade to. For example a script that upgrades to OSCAR version 0.13.0
should be named `0.13.0.jl`. 
There is also the possibility to have multiple upgrade scripts per version, this is to accommodate file serialized with DEV versions.
In this case the upgrades should be named `1.6.0-n.jl` where `n` is the `n`th upgrade in the sequence of upgrades that will upgrade a file to the `1.6.0` version.
To guarantee that upgrades occur in the correct order it is important that they are included (`include("/path/to/upgrade")`) in the correct order in `src/Serialization/Upgrades/main.jl`.

```@docs
Oscar.Serialization.UpgradeScript
```

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

# External Implementations

Any external body implementing a save/load following the `.mrdi` format specification
and using the OSCAR namespace should be sure to check validity against our schema defined
[here](https://www.oscar-system.org/schemas/mrdi.json).

We make no attempt whatsoever to verify the mathematics of the file, and neither
should anyone implementing a save/load. Loading should not throw a parse error
if the mathematics of the file is incorrect, the file should be parsed and allow
the computer algebra system to throw the error. We cannot guarantee that any file
that has been manipulated by hand is still valid and should be validated against
the schema. In the same way we cannot guarantee that any files created externally
are valid in terms of the mathematics either, these will not lead to a parse error
but instead will be handle as though the incorrect input has been passed to one
of the OSCAR functions.

External implementations should not be expected to read or write all possible OSCAR types.
It is perfectly valid for external implementations to throw parse errors when a certain
file format is unexpected. For example, OSCAR will parse a `QQFieldElem` that has data value
"0 0 7 // - 1 0" as `-7//10`, even though this is not how it is serialized. We feel
we should not restrict users when deserializing to formats that may have issues deserializing
the same format externally.

Allowing extensions to JSON is not recommended, this is to keep the scope
of possible software that can parse the given JSON as large as possible.
For example some JSON extensions allow comments in the files, OSCAR cannot
parse such JSONs and we recommend that any comments should be placed in the
meta field.

When writing UUIDs, adhere to version four UUIDs specified by RFC 4122.
