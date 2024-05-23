# [Serialization](@id dev_serialization)

We implement our serialization in accordance with the [MaRDI](https://www.mardi4nfdi.de/about/mission) file format specification described [here](https://arxiv.org/abs/2309.00465).
Which means we use a JSON extension to serialize data.

This document summarizes the serialization efforts of OSCAR, how it is supposed
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

# More information

   (1) We need name for the serialization language
      (e.g. I can say that CoCoALib has an "OpenMath" interface, but
       what do I say for Oscar-JSON-serialization?)

We serialize our data according to the MaRDI file format specification,
using the Oscar namespace. CoCoA would then serialize and deserialize
files in the mardi file format that use the oscar namespace.


(2)  We need to specify precisely what a receiver must accept:
      (2a) plain JSON syntax, or allow some extensions?
           [personally I'd like to allow end-of-line comments]

           I think we should be as strict as possible here to allow all
           possible json parsers to allow reading.
           We can put comments into the meta data section and still have a valid
           json format.
      (2b) structure of the "parse tree" (assuming that the byte stream
           was valid JSON syntax)

           is this not the schema?
           [here](https://www.oscar-system.org/schemas/mrdi.json)
           
      (2c) what behaviour is expected/permitted if the "parse tree" is
           not as expected?
          [probably "graceful degradation", i.e. error/exception;
           but should the receiver attempt to tell the sender if
           in an interactive context?]

       We currently don't validate against the schema, we could add a message about
       validating against the schema if there is an error on load.

(3)  Parse tree structure
      (3a)  We have a list of permitted "keys", and contexts in which they
            may appear -- this must be documented
      (3b)  Are unexpected keys permitted?  Or is that an error?

          this is fine, since they won't be accessed,
          
      (3c)  Are duplicate keys permitted?  If so, how to handle them?
            What if there are 2 or more "_ns" keys?

            Duplicate keys are not recomend in general JSON,
             [RFC 8259](https://www.rfc-editor.org/rfc/rfc8259#section-4)
             says they should be unique. 


      (3d)  Are the keys case-sensitive?

            yes 
      (3e)  May keys contain (leading/trailing) whitespace?
            no

(4)  UUIDs
      (4a)  Can these be any strings?  With max/min length?
            Or must they follow the variant 1 big-endian format
           (namely, hexadecimal with dashes)?
           no, version four UUIDs specified by RFC 4122
           
      (4b)  In "_refs" the value of each key has to be an object?
            In particular, we do not allow the value to be just another
           reference, otherwise an error must be signalled?

           This is handled by the schema. If it is not a an object with
           a _type key it is not valid
           
      (4c)  References must not create an infinite loop
            [receiver must detect this?]
            probably a good idea
            

(5)  Integer & rational literals
      (5a)  integer & rational literals are represented as non-empty strings
            [regardless of how many digits they have]
            yes
            
      (5b)  The only permited characters are an optional initial "-" and
            decimal digits; whitespace and other characters are not
permitted.
        Whitespace characters need to be escaped otherwise they are not valid json.
        
      (5c)  The integer zero is represented as "0";
            is "-0" also allowed?

            seems so, this is due to how julia parses the string "-0" to Int
            
      (5d)  Unnecessary leading zeroes are not allowed
            [is this too strict?]

            I (Antony) think we can maybe leave this up to the implementation?
            Why deny a user loading a file if it works in one software.
            The given implementation can spit out an error if it's unhappy
            about certain cases it doesn't feel like dealing with.
            Maybe we can be more precise by adding a warning about manually the files.
   
            
      (5e)  A rational literal is an integer literal optionally followed by
            a division-mark followed by a positive integer;
           we do NOT require that numerator and denominator be coprime
           [do we allow more than one division-mark?  Currently only "//"
            is permitted; JAA finds this weird, and thinks that "/" is
            more widely recognized]

            I (Antony) think the "//" makes sense since we are using Oscar and julia semantics.
            This is a hard a question and I think we should get more opinions on it.
            
      (5f)  Whitespace & other chars are not allowed in a rational literal.

            again I feel that anything that can be read into oscar shouldn't be stopped
            from loading into Oscar, .i.e we don't want to enforce how people write
            but instead recommend how they should, so as to remain interoperable
            using Oscar semantics with other systems that have their implementation.
            
      (5g)  Is zero denominator allowed?

            I don't think we should be checking the mathematics of the file. Currently
            we allow zero denominator in the file but this leads to a division error
            
      (5h)  Currently I'm definitely in favour of forbidding whitespace
            *inside* an integer/rational literal; this would mean that your
           implementation is not compliant!

           The point of the namespace is to say wtv can be read by oscar is valid.
           If your implementation wants to complain I think this also should be perfectly valid.
           
      (5j)  Do we allow a leading "+"?  What about "+0"?       

            +n seems to not be a valid big integer for any n
            seems to be ok for Int though. 

(6)  Other questions
      (6a)  the serialization of a matrix has type MatSpace(*,2,3) but the
            data array is of size 4-by-5.  What should happen?

            So in Oscar this will throw an error. But not a parse error,
            in terms of Oscar semantics the mathematics of the file is wrong.
            We should be clear in this distinction.
        
            
      (6b)  the serialization of a matrix has type MatSpace(*,4,5) but the
            data array is of size 2-by-3.  What should happen?

            same thing as above
            
      (6c)  the serialization of a matrix has non-rectangular data array;
            what should happen?

            parse error
            
      (6d)  the serialization of a prime finite field element has integer
            value outside the interval [-p+1, p-1].  Is this OK?

            yes
            
      (6e)  the serialization of a prime finite field element has a rational
            value whose denominator is coprime to p.  Is this OK?

            no
            
      (6f)  the serialization of an element in ZZ/nZZ has a rational
            value whose denominator is not coprime to n.  Is this OK?

            no

      (6g)  the serialization of a polynomial has negative exponent
            (or huge exponent); what should happen?
            exponents need to be positive ints
            
(7)  More sundries
  (a) May the "_refs" part contain unneeded bits?
      For instance, suppose CoCoA can handle only matrices of
integers/rationals,
      and indeed that is what the data is, BUT the "_refs" contains also a
      polynomial ring; must CoCoA accept this?  Or may it give an error?
      [here I assume that the data does not refer to bit about poly rings]

      implementations are free to intepret what they can and forget the rest,
      if cocoa is able to read a certain sub tree then it should not be stopped
      from doing so.

  (b) We discussed somewhat vaguely the notion of a "session"; if a session
      may involve the exchange of several messages then it may be useful to
      be able to specify new "_refs" even though some have already been
      specified.  Note that I am vague in my own mind here; so maybe what
      I have written is not very clear.  Perhaps the real question is:
      do we want to create a notion of "session" which comprises several
      JSON messages/objects ORDERED TEMPORALLY?

      no, sessions don't really exist, what happens is files from a session exist.
      loading any number of files from a specific session in any order should
      have no trouble loading and each of the respectively loaded objects should
      have no trouble with coercion provided none existed in the session in which
      they were saved in. Loading files from a session and saving them again should
      use the same uuids, it may happen that this is undesirable to the users
      some functionality for clearing the global serializer state of the session
      should be provided.

(8)

(1)  You have probably already thought of this...
      every JSON message must be either accepted by every receiver
      or rejected by every receiver (provided the receiver can
      "understand" the content).  Otherwise how do we achieve
      "reproducibility"?

      a good example here is some pdf viewers can load and play video,
      but it is not enforced that implementations should allow
      for playing video.

      same thing here

(1a) Here is an example of what we should not permit:
      An integer value is an optional sign followed by
      a string of decimal digits of length at least one;
      a receiver MAY CHOOSE TO allow leading/trailing
      whitespace. [["may choose to" is not allowed]]

      I think we can only say what is recommended and not try to
      strictly enforce struct, others we are creating too much overhead.
      Implementations are allowed to complain, but if it is able to
      load into Oscar then we should not prevent it. 
      
(2)  Either a "session" is a complete, indivisible JSON
      structure, or if it contains "separable" JSON structures
      then each of these structures must be individually complete
      (otherwise a separated JSON sub-stucture taken out of
      context may be complete enough to be comprehended); here
      "context" presumably means the "_refs"

      files are complete, unless they use a different serializer
      (we should also document the IPC Serializer?)
      Not sure if we really want to offer storing sessions,
      in which case storing sessions should just be a collection
      of files anyway.
      
(3)  Regarding matrices: there is an apparent redundancy because
      the matrix dimensions are given explicitly in the "MatSpace"
      type specification, and also given implicitly in the array
      of values for the matrix entries.  However giving the matrix
      dimensions in the "MatSpace" is necessary if the matrix has
      zero rows -- perhaps not necessary if there are zero columns.
      One could argue that if the matrix contains no entries then
      we can skip the "_data" subtree -- I'm unsure whether such an
      exceptional case is a good idea.

      not really sure what the question is here?

      we should allow saving and loading empty times as to not throw
      errors when people are trying to serialize their experiments
