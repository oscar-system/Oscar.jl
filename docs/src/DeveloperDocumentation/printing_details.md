# Details on printing in Oscar

The following dection contains more details and examples on how to implement
Oscar's 2+1 printing modes. The specifications and a minimal example may be
found in the [Developer Style Guide](@ref).


### Implementing show functions

Here is the translation between `:detail`, `one line` and `:supercompact`.

```
print(io, "text/plain", x)                 # detailed printing
print(io, x)                               # one line printing
print(IOContext(:supercompact => true), x) # supercompact printing
```

For reference, string interpolation `"$(x)"` will also use `print(io, x)`.

#### Mockup

##### Detailed printing with a new line.

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

##### Detailed printing in a single line.

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

##### The following is not working as expected and should not be used

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
