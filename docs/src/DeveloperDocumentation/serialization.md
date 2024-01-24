# Serialization

This document summarises the serialization efforts of OSCAR, how it is supposed
to work, how it works and the overall goal.
[Serialization](https://en.wikipedia.org/wiki/Serialization) broadly speaking
is the process of reading and writing data. There are many reasons for this
feature in OSCAR, but the main reason is communication on mathematics by
mathematicians.

## How it works
The mechanism for saving and loading is very simple. It is implemented via two
methods `save` and `load`, and works in the following manner:
```
julia> save("/tmp/fourtitwo.json", 42);

julia> load("/tmp/fourtitwo.json")
42

```
As hinted by the filename, OSCAR writes a file in JSON format. The file looks
as follows:
```
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
It contains the version of OSCAR used for serialization. The content is "42",
it represents a `Base.Int`, according to the `_type` field.


## Implementation
All files for serialization can be found in the folder `src/Serialization`. The
naming conventions of the files there follows the overall structure of OSCAR,
i.e. the file `src/Serialization/PolyhedralGeometry.jl` contains functions for
serializing objects of the polyhedral geometry section.

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
2. `save_object`, `load_object`
3. `save_type_params`, `load_type_params`

#### `save_type_object` / `load_type_object`

For the most part these functions should not be touched, they are high level
functions and are used to (de)serialize the object with its
type information as well as its data. The data and type nodes are
set in `save_typed_object` resulting in a "data branch" and "type branch".
The usage of these functions can be used inside `save_object` / `load_object`
and `save_type_params` / `load_type_params`. However using `save_typed_object` inside
a `save_object` implementation will lead to a verbose format and should at some
point be moved to `save_type_params`. Their implemention can be found in the
`main.jl` file.

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
```
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
```
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

```
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
```
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

```
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
```
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
```
function save_object(s::SerializerState, obj:NewType)
  save_object(s, obj.1)
end
```

While this will throw an error
```
function save_object(s::SerializerState, obj:NewType)
  save_object(s, obj.1, :key)
end
```

If you insist on having a key you should use a `save_data_dict`.
```
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

#### `save_type_params` / `load_type_params`

The serialization mechanism stores data in the format of a tree, with the
exception that some nodes may point to a shared reference. The "data branch"
is anything that is a child node of a data node, whereas the "type branch" is
any information that is stored in a node that is a child of a type node.
Avoiding type information inside the data branch will lead to a more
efficient serialization format. When the `uses_params` is set when
registering the type with [`@register_serialization_type`](@ref)
(de)serialization will use `save_type_params` / `load_type_params`
to format the type information.
In general we expect that implementing a `save_type_params` and
`load_type_params` should not always be necessary. Many types
will serialize their types in a similar fashion for example serialization
of a `FieldElem` will use the `save_type_params` / `load_type_params` from
`RingElem` since in both cases the only parameter needed for such types
is their parent.

### Serializers

The code for the different types of serializers and their states is found in the
`serializers.jl` file. Different serializers have different use cases, the
default serializer `JSONSerializer` is used for writting to a file. Currently
the only other serializer is the `IPCSerializer` which at the moment is
quite similar to the `JSONSerializer` except that it does not store the refs of
any types that are registered with the `uses_id` flag. When using the `IPCSerializer`
it is left up to the user to guarantee that any refs required by a process are sent
prior.

### Upgrades

All upgrade scripts can be found in the `src/Serialization/Upgrades` folder.
The mechanics of upgrading are found in the `main.jl` file where the
[`Oscar.upgrade`](@ref) function provides the core functionality. Upgrading
is triggered during [`load`](@ref) when the version of the file format
to be loaded is older than the current Oscar version.

```@docs
Oscar.upgrade
Oscar.upgrade_data
```

#### Upgrade Scripts

All upgrade scripts should be contained in a file named after the version
they upgrade to. For example a script that upgrades to Oscar version 0.13.0
should be named `0.13.0.jl`.

```@docs
Oscar.UpgradeScript
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
