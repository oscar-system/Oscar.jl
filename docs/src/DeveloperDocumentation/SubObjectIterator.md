# `SubObjectIterator`

Many of the objects in the field of *Polyhedral Geometry* mask a `BigObject`
from `Polymake.jl`. These big objects have properties which can easily be
accessed via julia's dot syntax. The return commonly does not adhere to the
mathematical or the typing conventions of `Oscar`; many properties encode
information about a collection of mathematical objects within a single data
object.

The `SubObjectIterator` is a precise and flexible tool to directly access and/or
process the desired properties of any `Polymake.BigObject`, but it requires
specific interface definitions to work properly for each context. The user can
thus profit from an easily understandable and usable iterator.

This guide is meant to communicate the application of the `SubObjectIterator`
for developers, utilizing existing code as reference and examples.

## Creating a working `SubObjectIterator`

The formal definition of the `SubObjectIterator` in `src/PolyhedralGeometry/iterators` is:

```julia
struct SubObjectIterator{T} <: AbstractVector{T}
    Obj::Polymake.BigObject
    Acc::Function
    n::Int
    options::NamedTuple
end
```

An instance can be created by passing values for all fields, while `options`
is optional.

Trivially, `Obj` is the `Polymake.BigObject` whose property is to be accessed.
The other fields will each be explained in an upcoming section.

### Length
As an `AbstractVector`, the `SubObjectIterator` has a length. Due to the nature
of `Polymake.BigObject`s this length is constant for any property. Sometimes the
length can easily be derived as a by-product of pre-computations when creating
an instance of `SubObjectIterator`. To avoid performing unnecessary computations
afterwards, the value is set at construction in `n`.

### Access function

Optimally retrieving and converting the elements varies strongly between the
contexts in which a `SubObjectIterator` is created. Thus its `getindex` method
redirects the call to the (internal) function `Acc`:

```julia
function Base.getindex(iter::SubObjectIterator{T}, i::Base.Integer) where T
    @boundscheck 1 <= i && i <= iter.n
    return iter.Acc(T, iter.Obj, i; iter.options...)
end
```

From this call we can see that the access function's signature needs to satisfy
certain requirements for the `SubObjectIterator` to work. The arguments are:

1. `T`: The return type.
2. `iter.Obj`: The `Polymake.BigObject` whose property is to be accessed.
3. `i`: The index.
4. `iter.options`: Additional arguments. Will be explained later.

Let us look at an example how we can utilize this interface. The following is
the implementation to access the rays of a `Cone`:

```julia
rays(as::Type{RayVector{T}}, C::Cone) where T = SubObjectIterator{as}(pm_object(C), _ray_cone, nrays(C))

_ray_cone(::Type{T}, C::Polymake.BigObject, i::Base.Integer) where T = T(C.RAYS[i, :])
```

Typing `r = rays(RayVector{Polymake.Rational}, C)` with a `Cone` `C` returns a
`SubObjectIterator` over `RayVector{Polymake.Rational}` elements of length
`nrays(C)` with access function `_ray_cone`. With the given method of this
function, `getindex(r, i)` returns a `RayVector{Polymake.Rational}` constructed from
the `i-th` row of the property `RAYS` of the `Polymake.BigObject`.

The user does never directly create a `SubObjectIterator`, so type restrictions
made where it is created can be assumed to hold. In our example `_ray_cone` will
always be called with `T<:RayVector`.

One can define several methods of the access function to ideally read and
process data. Consider `facets(as::Type{T}, C::Cone)`. Depending on the return
type we offer three methods:

```julia
_facet_cone(::Type{T}, C::Polymake.BigObject, i::Base.Integer) where T<:Union{Polyhedron, AffineHalfspace} = T(-C.FACETS[[i], :], 0)

_facet_cone(::Type{LinearHalfspace}, C::Polymake.BigObject, i::Base.Integer) = LinearHalfspace(-C.FACETS[[i], :])

_facet_cone(::Type{Cone}, C::Polymake.BigObject, i::Base.Integer) = cone_from_inequalities(-C.FACETS[[i], :])
```

#### Additional Methods

The `SubObjectIterator` can moreover be understood as a mathematical collection
the sense that one can

1. ask for specific information encoded in the data or
2. use this collection as an argument for construction another mathematical object.

The first case is covered by adding methods to specific internal functions.
Remember implementation of `rays` discussed above. It makes sense to define
a `vector_matrix` method on its output, encoding the rays of the cone as a
single matrix based on a convention applied throughout `Oscar`. The function's
implementation a user calls in this case is evaluated to these lines:

```julia
vector_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Rational}}) = matrix(QQ, Matrix{fmpq}(_vector_matrix(Val(iter.Acc), iter.Obj; iter.options...)))
vector_matrix(iter::SubObjectIterator{<:AbstractVector{Polymake.Integer}}) = matrix(ZZ, _vector_matrix(Val(iter.Acc), iter.Obj; iter.options...))
_vector_matrix(::Any, ::Polymake.BigObject) = throw(ArgumentError("Vector Matrix not defined in this context."))
```

Two functionalities are defined this way:

1. The call of `vector_matrix(iter)` is redirected to
   `_vector_matrix(Val(iter.Acc), iter.Obj)`. If that method is not defined
   for the value type of the access function, it falls back to throwing an
   error.
2. The matrix received from step 1 is converted from `Polymake.jl` format
   to `Oscar` format.

So by defining the following we have a fully functional `vector_matrix` method
in the context of `rays`:

```julia
_vector_matrix(::Val{_ray_cone}, C::Polymake.BigObject) = C.RAYS
```

The second case is solved with defining a special `_matrix_for_polymake` method.
One just hast to name the internal function that returns the desired matrix.
This way one has the ability to precisely control how the iterator works
internally in specific contexts, even if there happen to be multiple additional
matrix functions.

Again, the call `matrix_for_polymake(iter)` will either redirect to the defined
method or fall back to throwing an error if there is none:

```julia
function matrix_for_polymake(iter::SubObjectIterator)
    if hasmethod(_matrix_for_polymake, Tuple{Val{iter.Acc}})
        return _matrix_for_polymake(Val(iter.Acc))(Val(iter.Acc), iter.Obj; iter.options...)
    else
        throw(ArgumentError("Matrix for Polymake not defined in this context."))
    end
end
```

For `rays(C::Cone)` this reduces the implementation to the following line:

```julia
_matrix_for_polymake(::Val{_ray_cone}) = _vector_matrix
```

With `matrix_for_polymake` the output of `rays` can be handled as a usual matrix
and constructors or other functions can easily be extended by additionally
allowing `SubObjectIterator` as an argument type. E.g. the signature of one of
the `Cone` constructors now looks like this while the body has not changed:

```julia
Cone(R::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem, AbstractMatrix}, L::Union{SubObjectIterator{<:RayVector}, Oscar.MatElem, AbstractMatrix, Nothing} = nothing; non_redundant::Bool = false)
```

There also are `linear_matrix_for_polymake` and `affine_matrix_for_polymake`
used in the context of linear and affine halfspaces/hyperplanes. Defining this
functionality in a context works the same way as for `matrix_for_polymake`; you
can create a new method of `_linear_matrix_for_polymake` or
`_affine_matrix_for_polymake`. It suffices to define the most relevant of these
two; the other one will be derived, if possible. Also, `halfspace_matrix_pair`
is defined in terms of `affine_matrix_for_polymake`, so this does not need
another implementation.

The example code for `rays(C::Cone)` has covered every line of the
implementation by now, but we had different code in between, so let us summarize
and take a look at what the whole implementation actually looks like:

```julia
rays(as::Type{RayVector{T}}, C::Cone) where T = SubObjectIterator{as}(pm_object(C), _ray_cone, nrays(C))

_ray_cone(::Type{T}, C::Polymake.BigObject, i::Base.Integer) where T = T(C.RAYS[i, :])

_vector_matrix(::Val{_ray_cone}, C::Polymake.BigObject) = C.RAYS

_matrix_for_polymake(::Val{_ray_cone}) = _vector_matrix
```

### `options`

Sometimes you need further arguments to specify the returned data. These
arguments are set at construction of the `SubObjectIterator` and later passed to
the corresponding functions as keyword arguments.

A good example how to use this is `faces(C::Cone, face_dim::Int)`. It is not
enough to know that our `SubObjectIterator` is set in the context of faces of
cones; `face_dim` will be relevant for any type of access occurring in the future.

```julia
function faces(C::Cone, face_dim::Int)
   n = face_dim - length(lineality_space(C))
   n < 1 && return nothing
   return SubObjectIterator{Cone}(C.pm_cone, _face_cone, size(Polymake.polytope.faces_of_dim(pm_object(C), n), 1), (f_dim = n,))
end
```

When this method is called with meaningful input, it creates a
`SubObjectIterator` where the last argument is a `NamedTuple` specifying that
`f_dim = n`. The information encoded in this `NamedTuple` will be passed as
keyword arguments when calling the access function or any additional method
(reconsider their definitions). This allows us to directly ask for that data
when implementing these methods:

```julia
function _face_cone(::Type{Cone}, C::Polymake.BigObject, i::Base.Integer; f_dim::Int = 0)
   return Cone(Polymake.polytope.Cone(RAYS = C.RAYS[collect(Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(C, f_dim)[i])), :], LINEALITY_SPACE = C.LINEALITY_SPACE))
end

function _ray_indices(::Val{_face_cone}, C::Polymake.BigObject; f_dim::Int = 0)
   f = Polymake.to_one_based_indexing(Polymake.polytope.faces_of_dim(C, f_dim))
   return IncidenceMatrix([collect(f[i]) for i in 1:length(f)])
end
```

## Extending the interface

The additional methods offer an intuitive way of interaction for the user, but
their current selection is not carved in stone. You can easily add more similar
methods by extending the list that is iterated over to generate the code. Which
list that is usually depends on the output format. `vector_matrix` returns
matrices with either integer or rational elements. The same capabilities hold
for `point_matrix` and `generator_matrix`:

```julia
for (sym, name) in (("point_matrix", "Point Matrix"), ("vector_matrix", "Vector Matrix"), ("generator_matrix", "Generator Matrix"))
    M = Symbol(sym)
    _M = Symbol(string("_", sym))
    @eval begin
        $M(iter::SubObjectIterator{<:AbstractVector{Polymake.Rational}}) = matrix(QQ, Matrix{fmpq}($_M(Val(iter.Acc), iter.Obj; iter.options...)))
        $M(iter::SubObjectIterator{<:AbstractVector{Polymake.Integer}}) = matrix(ZZ, $_M(Val(iter.Acc), iter.Obj; iter.options...))
        $_M(::Any, ::Polymake.BigObject) = throw(ArgumentError(string($name, " not defined in this context.")))
    end
end
```

The second string (`name`) of each pair determines the name that is printed in error
messages.

If required, one can of course write completely new functions to extend the
interface.
