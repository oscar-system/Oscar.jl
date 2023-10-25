# Details on printing in Oscar

The following dection contains more details and examples on how to implement
Oscar's 2+1 printing modes. The specifications and a minimal example may be
found in the [Developer Style Guide](@ref).


## Implementing show functions

Here is the translation between `:detail`, `one line` and `:supercompact`,
where `io` is an `IO` object (such as `stdout` or an `IOBuffer`):

```
show(io, MIME"text/plain"(), x)                # detailed printing
print(io, x)                                   # one line printing
print(IOContext(io, :supercompact => true), x) # supercompact printing
```

For reference, string interpolation `"$(x)"` uses one line printing via `print(io, x)`,
while on the REPL detailed printing is used to show top level objects.

!!! warning "display"
    Please do not use `display`! From the [Julia documentation of
    `display`](https://docs.julialang.org/en/v1/base/io-network/#Base.Multimedia.display):
    "In general, you cannot assume that `display` output goes to `stdout`
    [...]". In particular, the output of `display` will not work in the
    `jldoctest`s.

### Mockup

#### Detailed printing with a new line

```julia
struct NewRing
  base_ring
end

base_ring(R::NewRing) = R.base_ring
```
The following is a template for `detailed` printing.
Note that at least one new line is needed for technical reasons.
see below why.

```julia
function Base.show(io::IO, ::MIME"text/plain", R::NewRing)
  println(io, "I am a new ring")  # at least one new line is needed
  println(io, "I print with newlines")
  print(io, base_ring(R)) # the last print statement must not add a new line
end
```

The following is a template for `one line` and `:supercompact` printing.
```julia
function Base.show(io::IO, R::NewRing)
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "supercompact printing of newring ")
  else
    # nested printing allowed, preferably supercompact
    print(io, "one line printing of newring with ")
    print(IOContext(io, :supercompact => true), "supercompact ", base_ring(R))
  end
end
```
And this is how it looks like:
```julia
julia> R = NewRing(QQ)
I am a new ring
I print with newlines
QQ

julia> [R,R]
2-element Vector{NewRing}:
 one line printing of newring with supercompact QQ
 one line printing of newring with supercompact QQ

```

#### Detailed printing in a single line

This version needs to be used in case the detailed
printing does not contain newlines.
Then detailed and one line printing agree.
The `if` clause takes care of supercompact printing as well.

```julia
struct NewRing2
  base_ring
end

base_ring(R::NewRing2) = R.base_ring

function Base.show(io::IO, R::NewRing2)
  if get(io, :supercompact, false)
    # no nested printing
    print(io, "supercompact printing of newring")
  else
    # nested printing allowed, preferably supercompact
    print(io, "I am a new ring and always print in one line " )
    print(IOContext(io, :supercompact => true), base_ring(R))
  end
end
```
And this is how it looks like:
```julia
julia> R = NewRing2(QQ)
I am a new ring and always print in one line QQ

julia> [R,R]
2-element Vector{NewRing2}:
 I am a new ring and always print in one line Rational Field
 I am a new ring and always print in one line Rational Field

julia> print(IOContext(Base.stdout, :supercompact => true) ,R)
supercompact printing of newring
```

The `supercompact` printing uses an `IOContext` (see [IOContext](https://docs.julialang.org/en/v1/base/io-network/#Base.IOContext) from
the Julia documentation) to pass information to other `show` methods invoked
recursively (for example in nested printings). The same mechanism can be used to
pass other context data. For instance, this is used by the `Scheme` code in
some nested printings which invoke several objects whose printing depends on a
given *covering*: we use `IOContext` to pass a fix covering to the printing
of each sub-object for consistency and readability.

#### The following is not working as expected and should not be used

This example does not work correctly because the `detailed` printing does not
include a newline, which is expected by the Julia printing system. To correctly
support single line `detailed` printing, read the preceding section.

```julia
function Base.show(io::IO, ::MIME"text/plain", R::NewRing)  # do not implement me like this
  print(io, "I am a new ring with a detailed printing of one line")
end
```

Then the following will *not* be used for array/tuple printing.
It will be used for `print(io, R::NewRing)` though.

```julia
function Base.show(io::IO, R::NewRing)
  if get(io, :supercompact, false)
    print(io, "supercompact printing of newring")
  else # this is what we call one line
    print(io, "one line printing of newring with ")
    print(IOContext(io, :supercompact => true), "supercompact ", R.base_ring)
  end
end
```
This example illustrates the unexpected behavior.
```julia
julia> R = NewRing(1)

julia> R
I am a new ring with a detailed printing of one line

julia> [R,R]  # one line printing is ignored
2-element Vector{NewRing}:
 I am a new ring with a detailed printing of one line
 I am a new ring with a detailed printing of one line

julia> print(Base.stdout, R)
one line printing of newring with supercompact QQ
```

## Advanced printing functionality

To facilitate printing of nested mathematical structures, we provide a modified
`IOCustom` object. To create one, we use the following command:

```@docs
AbstractAlgebra.pretty(::IO)
```

The `IOCustom` object allows one to locally control:
- indentation using `Indent()` and `Dedent()`,
- capitalization using `Lowercase()` and `LowercaseOff()`.

### Example

We illustrate this with an example

```
struct A{T}
  x::T
end

function Base.show(io::IO, a::A)
  io = AbstractAlgebra.pretty(io)
  println(io, "Something of type A")
  print(io, AbstractAlgebra.Indent(), "over ", AbstractAlgebra.Lowercase(), a.x)
  print(io, AbstractAlgebra.Dedent()) # don't forget to undo the indentation!
end

struct B
end

function Base.show(io::IO, b::B)
  io = AbstractAlgebra.pretty(io)
  print(io, AbstractAlgebra.LowercaseOff(), "Hilbert thing")
end
```

At the REPL, this will then be printed as follows:
```
julia> A(2)
Something of type A
  over 2

julia> A(A(2))
Something of type A
  over something of type A
    over 2

julia> A(B())
Something of type A
  over Hilbert thing
```

Moreover, one can control the pluralization of nouns when printing a set of
elements with a variable number of objects. For this, one can use `ItemQuantity`:

### Example

We illustrate this with an example

```
julia> struct C{T}
       x::Vector{T}
       end

julia> function Base.show(io::IO, c::C{T}) where T
       x = c.x
       n = length(x)
       print(io, "Something with ", AbstractAlgebra.ItemQuantity(n, "element"), " of type $T")
       end
```

At the REPL, this will then be printed as follows:
```
julia> C(Int[2,3,4])
Something with 3 elements of type Int64

julia> C(Int[])
Something with 0 elements of type Int64

julia> C(Int[6])
Something with 1 element of type Int64
```

## LaTeX and Unicode printing

### LaTeX output
Some types support LaTeX output.
```
julia> Qx, x = QQ["x"];

julia> show(stdout, "text/latex", x^2 + 2x + x^10)
x^{10} + x^{2} + 2 x

julia> show(stdout, "text/latex", Qx[x x^2; 1 1])
\begin{array}{cc}
x & x^{2} \\
1 & 1
\end{array}
```

```
Base.show(io::IOContext, ::MIME"text/latex")
```

### Unicode printing
Per default output should be ASCII only (no Unicode). Implementors of
`Base.show` and related functions can branch on the output of
`Oscar.is_unicode_allowed()` to display objects using non-ASCII characters.
This will then be used for users which enabled Unicode using
`allow_unicode(true)`. Note that

- there must be a default ASCII only output, since this is the default setting
  for new users, and
- OSCAR library code is not allowed to call `Oscar.allow_unicode`.

Here is an example with and without output using Unicode:

```julia
  struct AtoB
  end

  function Base.show(io::IO, ::AtoB)
    if Oscar.is_unicode_allowed()
      print(io, "A→B")
    else
      print(io, "A->B")
    end
  end
```
