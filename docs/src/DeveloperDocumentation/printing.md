# Printing in Oscar

## The 2+1 new print modes
We propose two user print modes `:details` and `:oneline` and one internal print mode `:supercompact`. The latter is for use during recursion, e.g. to print the `base_ring(X)` when in `:oneline` mode. It exists to make sure that `:oneline` stays compact and human readable.

Top-level REPL printing of an object will use `:details` mode by default
 ```julia
 julia> X
 :details
 ```
 Inside nested structures, e.g. inside a `Vector`, the `:oneline` mode is used.
 ```julia
 julia> [X,X]
 3-element Vector{TypeofX{T}}
 :oneline
 :oneline
 :oneline
 ```

- `:details`; no restrictions except:
  - the output must make sense as a standalone without context to non-specialists
  - the number of output lines must fit in the terminal
  - the row length should fit the terminal (as it is for julia matrices);
  - if the object is simple enough use only one line
- `:oneline`;  the output must print in one line
  - should make sense as a standalone without context
  - Variable names/generators/relations should not be printed only their number.
- `:supercompact`; no nested printing
  - a single term or a symbol/letter mimicking mathematical notation
  - should usually only depend on the type and not of the type parameters or of the concrete instance - exceptions of this rule are possible e.g. for `GF(2)`
  -  no nested printing. In particular variable names and `base_ring` must not be displayed. This ensures that `:oneline` and `:supercompact` stay compact even for complicated things. If you want nested printing use `:oneline`. See it as a user readable version of the main (mathematical) type.

Every type must implement `:oneline` printing, and may optionally implement `:details` or `:supercompact`. The fallback is always `:oneline`.

One-line printing is used for substructures, e.g. the base ring, domain/codomain or in arrays, being in one line is important for formatting.

Finally `:supercompact` is to be used for information coming from substructures or subsubstructures e.g. the `base_ring` of the `base_ring` will be printed `:supercompactly.`


## Examples
For each example we fix a single parent and display its three print modes consecutively each in its own block
`:details`
`:oneline`
`:supercompact`
### Coefficient rings
```
Integer ring
```
```
Integer ring
```

```
ZZ
```
-----
```
Rational field
```
```
Rational field
```
```
QQ
```
----
##### Prime field
```
Galois field with characteristic 2
```
```
Galois field with characteristic 2
```
```
GF(2)
```
##### Finite field
```
Galois field of degree 7 over GF(29)
```
```
Galois field of degree 7 over GF(29)
```
```
GF(29^7)
```
##### Integer mod rings
```
Integers modulo 125
```
```
Integers modulo 125
```
```
ZZ/(125)  # ?
```
#### Numberfields
#### A simple number field
```
Number field of degree 10
with defining polynomial x^10 - 3x^9 + 7x^8 - 11x^7 + 13x^6 - 12x^5 + 9x^4 - 5x^3 + 3x^2 - 2x + 1
  over QQ
```
```
Number field of degree 10 over QQ
```
```
Number field
```
##### Fancy number field
```
Number field of degree 10
with defining polynomial x^10 - 3x^9 + 7x^8 - 11x^7 + 13x^6 - 12x^5 + 9x^4 - 5x^3 + 3x^2 - 2x + 1
  over
  Number field of degree 4
  with defining polynomial 9x^4 - 5x^3 + 3x^2 - 3x + 1
    over
    Number field of degree 2
    with defining polynomial x^2+1
    over
      QQ
```
```
Number field of degree 2 over Number field
```
```
Number field
```
----
### Groups

##### Permutation group
```
Permutation group of degree 7
```

```
Permutation group of degree 7
```

```
Permutation group
```
##### Polycyclic group
```
Polycyclic group of infinite order
```

```
Polycyclic group of infinite order
```

```
Polycyclic group
```
##### General linear group over finite field
```
General linear group of degree 24
  over Galois field of degree 7 over GF(29)
```
```
General linear group of degree 24 over GF(29^7)
```
```
General linear group
```

##### General linear group over number field
```
General linear group of degree 24
  over Number field of degree 2 over Number field
```
(note the oneline but not supercompact representation of the number field)
```
General linear group of degree 24 over Number field
```
```
General linear group
```
#### Matrix groups
```
Marix group of degree 24 over Rational field
```
```
Matrix group of degree 24 over QQ
```
```
Matrix group
```

### Rings
#### Common polynomial rings
```julia
Multivariate polynomial ring in x, y, z over QQ
Multivariate polynomial ring in 30 variables x1, x2, x3, ... , x30 over QQ
```
```julia
Multivariate polynomial ring in 30 variables over QQ
```

```
Multivariate polynomial ring
```
###### Examples for Polynomial rings
```julia
julia> P,x = polynomial_ring(QQ,[:x,:y,:z])
(Multivariate polynomial ring in 3 variables, QQMPolyRingElem[x, y, z])

```

#### Generic polynomial rings
```
Multivariate polynomial ring
  over Fraction field of Univariate polynomial ring over QQ
in 30 variables x1, x2, x3, x4, x5, x6, ..., x30
```

```
Multivariate polynomial ring in 30 variables over Univariate function field
```
```
Multivariate polynomial ring
```

#### Group Rings
Note the `:supercompact` printing in `:oneline` printing
```julia
Group ring
  of Polycyclic group of infinite order
  over Number field of degree 4
```
```julia
  Group ring of Polycyclic group over Number field
```
```julia
Group ring
```


---
#### Ideals
##### Polynomial ring
```
Ideal
  of Multivariate polynomial ring in 30 variables over QQ
generated by
 x1
 x2
 x5^23
   ⋮
 x30 + x2 +  some very big polynomial fitting several pages but which has no line break
```

```
Ideal with 126 generators in Multivariate polynomial ring
```

```
Ideal
```
#### Maximal order of a number field
I am not quite sure how to deal with the fact that there are always 2 generators and the first one is an integer. One could put them in one line only, but it would be "against the rule"
```julia
Ideal
  of Maximal order of Cyclotomic field of order 19
norm: 1123272193876636254792988373
minimum: 1013
two normal wrt: 1013
generated by
 1013
 z_19^9 + 389*z_19^8 + 1011*z_19^7 + 626*z_19^6 + 391*z_19^5 + 386*z_19^4 + 623*z_19^3 + 2*z_19^2 + 388*z_19 + 1012
```
Clearly, the number of generators is useless in this context. It is always 2. What would be useful? It seems to me like $I \cap \mathbb{Z}$?
```julia
julia> F,z = cyclotomic_field(19);

julia> OF = maximal_order(F);

julia> support(30*OF)
2-element Vector{NfOrdIdl}
 Ideal in maximal order
 Ideal in maximal order
 Ideal in maximal order
```

```julia
Ideal
```

#### Quotient Rings

#### A common example
```
Quotient
  of Multivariate polynomial ring in 7 variables over GF(2)
  by Ideal with 126 generators
```

```
Quotient of multivariate polynomial ring in 3 variables over GF(2)
```

```
Multivariate quotient ring
```
#### A fancy example

```
Quotient
  of Multivariate polynomial ring in 6 variables over Univariate function field
  by Ideal with 126 generators

```

```
Quotient of multivariate polynomal ring in 15 variables over Univariate function field
```

```
Multivariate quotient ring
```
### Spectra
```
Spectrum of
  Quotient of Multivariate polynomial ring in 5 variables over GF(2)
  by Ideal with 126 generators
```
```
Spectrum of Multivariate quotient ring with 3 relations over GF(2)
```

```
Spec of Multivariate quotient ring
```

### Maps
#### Maps of Rings
```
Ring homomorphism
  from Multivariate polynomial ring in 30 variables over QQ
  to Quotient of Multivariate polynomial ring in 3 variables over QQ.
defined by
 x1 -> x
 x2 -> y
 x3 -> z
 x4 -> 0
   ⋮
x30 -> 0
```

```
Hom: Multivariate polynomial ring -> Multivariate quotient ring
```

```
Ring homomorphism
```
#### Maps of Groups
```
Group homomorphism
  from Polycyclic group of size 24 with 3 generators
  to Matrix group of degree 24 with 10 generators over QQ
```
```
Hom: Polycyclic group -> Matrix group
```
```
Group homomorphism
```
#### Maps of rings but not of algebras
```julia
Ring homomorphism
  from Multivariate polynomial ring in 2 variables over Quadratic Field
  to  Multivariate polynomial ring in 3 variables over Quadratic field
  with coefficient map Automorphism of quadratic Field
defined by
 x -> x1
 y -> x1*x2
```

```julia
Hom: Multivariate polynomial ring -> Multivariate polynomial ring
```

```julia
Ring homomorphism
```

## Examples of Interactive use
### Primary decomposition
```julia
julia> p, (x,y,z) = polynomial_ring(GF(2), [:x, :y, :z])
(Multivariate Polynomial Ring in 3 variables over GF(2), fpMPolyRingElem[x, y, z])

julia> I = ideal([x^2 + y^2, z^3 - 1])
Ideal of Multivariate Polynomial Ring in 3 variables over GF(2)
generated by
 x^2 + y^2
 z^3 - 1

julia> PD = primary_decomposition(I)
2-element Vector{Tuple{MPolyIdeal{fpMPolyRingElem}, MPolyIdeal{fpMPolyRingElem}}}:
 (Ideal with 2 generators, Ideal with 2 generators)
 (Ideal with 2 generators, Ideal with 2 generators)

julia> PD[1][1]
Ideal of Multivariate Polynomial Ring in 3 variables over GF(2)
generated by
 z + 1
 x^2 + y^2
```

## Styleguide
### :details
```
Mathematical Type (aka supercompact)
  of $compact(Subobject1)
  in/over $oneline(subobject2)
generated by
$(some array)
```
### :oneline
- Only the first word is capitalized e.g. `Polynomial ring`
- one should just use `:supercompact` for nested printing in compact
```
Structure in $(ngens(X)) generator$(ngens(X)==1? "":"s") contained in $:supercompact(ambient(X)) over $:supercompact(base_ring(X))
```
- nested calls to `:oneline` (if you think them really necessary) must be at the end,
  so that one can read sequentially. Calls to `:supercompact` can be anywhere e.g.

  ```
  Structure in $(ngens(X)) generator$(ngens(X)==1? "":"s") $:supercompact(ambient(X)) over $:oneline(base_ring(X))
  ```
  because this
  ```
  Structure in 3 generators with 25 relations over  bla in 6 generators over blub in 1 generator
  ```
  is easier to read than
  ```
  Structure over bla over blub in 5 generators in 6 generators in 7 generators with 25 relations
  ```





# Questions, discussion, TODOs

## Grammar
In GAP we recently spent quite some effort to get away from printing things like
`in 1 variables` changing it to `in 1 variable`.
This might seem minor, and programmers like to ignore it, but it sounds quite jarring.
Dealing with this is not rocket science, we have a nice helper in GAP
that semi-automates this for us in 99% of all cases.

### Univariate vs Multivariate:
Since univariate and multivariate rings behave very differently, we want to distinguish them in supercompact printing:

`Univariate polynomial ring`
`Univariate quotient ring`
`Multivariate polynomial ring`
`Graded multivariate polynomial ring`
`Quotient of multivariate polynomial ring`
`Quotient of graded multivariate polynomial ring` ?


## General comments


## Unicode
Some examples for oneline printing using unicode
```julia
Polynomial ring in 10 variables over ℚ

Polynomial ring in 10 variables over 𝔽₂

Polynomial ring in 10 variables over 𝔽₃₁
```
The spacing of super and subscripts seems to depend on the font.
Here it looks great! But in my terminal it looks really bad.
But setting my font to monospace regular it looks really nice, just like here.
```julia
The hyperplane of ℝ² described by x₂ = 0
```

# Technical aspects

## How julia prints things

There are lots of open issues regarding with the printing in julia. I am too lazy to find those.

The important "features" for us are the following:
- `:oneline` is never set when elements of arrays or vectors are printed. Only for multi-dimensional arrays.
- The normal `X]enter` command in the REPL will always call `show(::IO, ::MIME"text/plain", X)` (which has some fallbacks, but this is the method we should overwrite for our "detailed" version)
- The printing of vector or Tuple elements is handled as follows:
    - If the `::MIME"text/plain"` `show` method contains a newline, it will use the `show(::IO, X)` method.
    - If there is no newline, use `show(::IO, ::MIME"text/plain", X)`.

  (This is of course not documented, but see the example below)
- So here is the upshot. Inside a `show(...)` call, one can never know if the caller is an array/vector/tuple printing method. If our `:detail` version contains no newline, the `:detail` version will also be used in array/tuple printing and the `:oneline` version will never be invoked.

Note that the symbols `:detail` or `:oneline` are just used in this document to differentiate between the different printing versions. They do not appear in any implementation, since this is not how printing works in julia.

## Summary

Here is the translation beetween `:detail`, `:oneline` and `:supercompact`.

```
print(io, "text/plain", x)                 # detailed printing
print(io, x)                               # oneline printing
print(IOContext(:supercompact => true), x) # supercompact printing
```

For reference, string interpolation `"$(x)"` will also use
```
print(io, x)
```

## Mockup

```julia
struct NewRing
  base_ring
end
```

```julia
function Base.show(io::IO, ::MIME"text/plain", R::NewRing)
    println(io, "I am a new ring")
    println(io, "I print with newlines")
    print(io, R.base_ring)
end
```

```julia
function Base.show(io::IO, R::NewRing)
    if get(io, :supercompact, false)
        print(io, "supercompact printing of newring ")
        print(io, R.base_ring) # this prints the :oneline version of newring
    else
        print(io, "oneline printing of newring with ")
        print(IOContext(io, :supercompact => true), "supercompact ", R.base_ring)
    end
end
```


### Alternative

Here is an alternative version, which needs to be used in case the detailed printing does not contain newlines.

```julia
function Base.show(io::IO, R::NewRing)
    if get(io, :supercompact, false)
      print(io, "supercompact printing of newring")
    else
      println(io, "I am a new ring and always print in one line: )
      print(IOContext(io, :supercompact => true), base_ring(R))
    end
end
```

### Not working as expected

```julia
function Base.show(io::IO, ::MIME"text/plain", R::NewRing)
    print(io, "I am a new ring with a detailed printing of one line")
end
```

Then the following will *not* be used for array/tuple printing. It will be used for `print(io, R::NewRing)` though.

```julia
function Base.show(io::IO, R::NewRing)
    if get(io, :supercompact, false)
        print(io, "supercompact printing of newring")
    else # this is what we call :oneline, but wrong terminology for julia
        print(io, "compact printing of newring with ")
        print(IOContext(io, :supercompact => true), "supercompact ", R.base_ring)
    end
end
```

## Some fallback functions

```julia
show_details(x) = show(stdout, "text/plain", x)

show_one_line(x) = show(stdout, x)

show_supercompact(x) = show(IOContext(stdout, :supercompact => true), x)
```

and also provide some macros

```
@show_details x
@show_one_line x
@show_super_compact x
```

## How we know all this about the printing system

```julia
julia> struct A end

julia> Base.show(io::IO, ::MIME"text/plain", a::A) = print(io, get(io, :compact, false), " A with text/display normal")

julia> A()
false A with text/display normal

julia> [A()]
1-element Vector{A}:
 false A with text/display normal

julia> Base.show(io::IO, ::MIME"text/plain", a::A) = print(io, get(io, :compact, false), " A with\ntext/display normal")

julia> [A()]
1-element Vector{A}:
 A()

julia> Base.show(io::IO, a::A) = print(io, get(io, :compact, false), " A normal")

julia> [A()]
1-element Vector{A}:
 false A normal

julia> (A(), A())
(false A normal, false A normal)

julia> [A() A(); A() A()]
2×2 Matrix{A}:
 true A normal  true A normal
 true A normal  true A normal

 julia> print(stdout, A())
false A normal

julia> "$(A())"
"false A normal"
 ```


 ### TODOs
 - Helper for `with 1 generator` and not `with 1 generators`
 - make indentation work by using stuff from https://github.com/KristofferC/IOIndents.jl
   the nested call will then be `print(indented(io), R.substruct)`
 - Max wants to add a few more examples so we can discuss them:
   - detailed group homomorphisms given by generator-images, including the mappings (`f1 -> (1,2)(3,4,5)`)
   - permutations
   - matrices (in detailed view they could print their ring?)
   - matrices-as-group-elements (should be distinguishable from plain matrices? at least in detailed mode)
 - Max: Right now all output starts with an uppercase letter. Drawback: this makes it harder to compose outputs.
   Perhaps we need an option to control that? Or we write a helper function which changes the case of the first letter of a string --
   but then we need to first print into strings, instead of directly printing to the output io, which means overhead and has other concerns
