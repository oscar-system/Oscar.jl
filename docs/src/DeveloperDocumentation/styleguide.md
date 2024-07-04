# Developer Style Guide

In general we aim to follow the [Julia Style Guide](https://docs.julialang.org/en/v1/manual/style-guide/)
but there are some exceptions due to our specific needs and a different background.

The content of this page are merely *guidelines*. There may be good reasons to
deviate from them in some cases; in that case just do so.

## General styleguide

- Use Julia conventions where applicable and when they don't contradict our
  own rules above.
- Unless really really necessary, don't add new dependencies. Every new
  dependency complicates the development workflow, in that we will need to stay
  compatible with this package. 
- If already existing types in OSCAR are almost what you need, consider
  improving them instead of writing your own. While it might be tempting to
  create a new polynomial ring type for the new application because some
  feature is missing, it causes a lot of work and compatibility issues: Will
  the new type support
  - normal functions (gcd, factor),
  - quotient fields,
  - modules and residue rings,
  - conversion to and from other already existing types?
- Whenever functions return the same mathematical object, but in different
  mathematical categories, the first argument should be the desired return
  type. One example is `projective_space(NormalToricVariety, *)` vs
  `projective_space(ProjectiveScheme, *)`. However, if the return type is
  different, even if the result describes the same mathematical object, it
  should be indicated in the function name, for example `automorphism_group` vs
  `automorphism_group_generators` vs `automorphism_list`.
- Whenever functions expect a ring, field, algebra, etc. as input they should
  be passed as the first argument, for example, `polynomial_ring(QQ, "x")`.
- Follow the mathematics. If your function needs a list of points, you should
  create a point-type (or use the one already there) and then use this.
  For user-facing functions, please do not use re-purposed lists, arrays,
  matrices...
- Input sanity checks should be enabled by default, they can then be disabled
  internally if they are known to be true, and manually by users.


## Naming conventions

The usual [Julia naming conventions](https://docs.julialang.org/en/v1/manual/style-guide/#Use-naming-conventions-consistent-with-Julia-base/)
apply to OSCAR, too (that said, for various reasons our code still violates
quite some of them; but in general we strive to reduce these).
Here is a summary of the naming convention followed in OSCAR:

- Use `CamelCase` for types and `snake_case` for *everything* else. (Internal
  functions do not have to follow these rules.) Types (and their constructor)
  tend to be in `CamelCase`. However, please *also* provide the constructor (or a
  constructor) in `snake_case`. As a user one usually does not know if something
  is a constructor or a function.
- For filenames we recommend using `snake_case.jl`.
- Noteworthy difference to Julia base is that we do not have exceptions for
  `is*` or `has*` as prefix.
  It is `is_foo` instead of `isfoo` and `has_bar` instead of `hasbar`.
  The main reason is to avoid awkward constructions like `isvery_ample`, while
  also being consistent.
  For compatibility with standard Julia, while staying consistent internally,
  we also provide aliases (using `AbstractAlgebra.@alias`) for various standard
  Julia functions, e.g. `is_one` as alias for `isone`
- For generic concepts choose generic names, based on general algebraic
  concepts, preferably not special names from your area of speciality.
- **Avoid direct access to members of our objects.** This means, do not use
  something like `A.foo`, instead use a suitable getter `get_foo(A)`, and if
  there is none, please write one or request that one be written. Internal
  member names are free to change at any time, but functions can be deprecated
  properly.
- In Julia we have multiple dispatch, so we do not need functions like
  `point_from_matrix` as the "from" part is clear by the type of the argument.
  It should be called `points(T::Matrix)` in some variation.
  Similarly for `matrix_to_points`. Of course it is fine to use them
  internally, where useful.


## Code formatting

### Editor configuration

Please check if your editor can be configured to honor our `.editorconfig`
file, see <https://editorconfig.org> for more information about this.

### JuliaFormatter

There is a `.JuliaFormatter.toml` in our git repository. To format your files,
first add `JuliaFormatter.jl` in Julia and then use
```
using JuliaFormatter
format_file("path/to/file/file.jl")
```

### Unicode

As most modern programming languages, Julia allows the use of Unicode, e.g.,
`α`, in the REPL as well as in source code. As this reduces accessibility
to various groups of users and developers, the use of Unicode should be kept
to a minimum. Here is a general principle:

> Do not use Unicode characters inside functions. See
> [Unicode printing](@ref) for the exception concerning printing.

### Whitespace

- Do not use tabs.
- Do not put spaces "inside" parenthesis.
- Do put spaces after commas.

Good example:
```julia
f(x, y) = x + 1
print(f(1, 2))
```

Bad example:
```julia
f( x,y ) = x + 1
print( f ( 1,2 ) )
```

### Loops and other control structures

- `for` loops should use `in` not `=`
- don't put spaces around the `:` in a range

Good example:
```julia
for i in 1:3
  println(i)
end
```

Bad example:
```julia
for i = 1 : 3
  println(i)
end
```

## Code structure

- do not nest loops and `if` clauses too deeply; if you are using 5 or more
  levels, then in general that's a hint that you should refactor; e.g.
  - by moving parts of the code into a separate function
  - by replacing guard constructs like
    ```julia
    for i in A
      if flag
        ...
      end
    end
    ```
    by
    ```julia
    for i in A
      if !flag
        continue
      end
      ...
    end
    ```
    or
    ```julia
    for i in A
      flag ||continue
      ...
    end
    ```
  - by merging loops: you can replace
    ```julia
    for i in A
      for j in B
        ...
      end
    end
    ```
    by
    ```julia
    for i in A, j in B
      ...
    end
    ```

- Functions should not have too many arguments.
  If you need a bunch arguments, chances are that introducing a new type
  makes it more readable.

- Functions should not be too long; very long functions are in general harder
  to understand; it is also more difficult to see all the code at once. Consider
  splitting the function into multiple ones, if it is sensibly possible.

- Every export statement must be confined to a single line; the intention is to
  make it easy to use tools like `git grep` to find exports. In general it is
  recommended export exactly one identifier per export statement. Exceptions may
  be made for certain tightly related identifiers, e.g. `is_finite`, `set_is_finite`
  and `has_is_finite` could be put on a single line. In general if multiple
  export statements appear in sequence, they must be sorted alphabetically.

However, as always, rules sometimes should be broken.

## Optional arguments for parents of return values

Several objects in OSCAR have `parent`s, e.g. polynomials, group elements, ... 
Whenever a function creates such objects from an input which does not involve 
the output's parent, we strongly recommend that 
the user should have the possibility to pass on this 
parent as a keyword argument under the name `parent`.
Beyond that you can make more entry points for such parents available for 
the user's convenience.

Let's see an example. Say, you want to implement the characteristic 
polynomial of a matrix. You could do it as follows:
```julia
  function characteristic_polynomial(A::MatrixElem)
    kk = base_ring(A)
    P, x = kk["x"]
    AP = change_base_ring(P, A)
    return det(AP - x*one(AP))
  end
```
You can see that the polynomial ring `P`, i.e. the parent of the output, 
is newly created in the body of the function. In particular, calling this 
function two times on two different matrices `A` and `B` might produce 
incompatible polynomials `p = det(A - x*one(A))` and `q = det(B - x*one(B))` 
with different parents. Calling `p + q` will result in an error. 

To solve this, we should have implemented the function differently:
```julia
  # Implementation of the recommended keyword argument signature:
  function characteristic_polynomial(
      A::MatrixElem;
      parent::AbstractAlgebra.Ring=polynomial_ring(base_ring(A), :t)[1]
    )
    AP = change_base_ring(parent, A)
    x = first(gens(ring))
    return det(AP - x*one(AP))
  end

  # Optional second signature to also allow for the specification of the 
  # output's parent as the first argument:
  function characteristic_polynomial(
      P::PolyRing,
      A::MatrixElem
    )
    coefficient_ring(P) === base_ring(A) || error("coefficient rings incompatible")
    return characteristic_polynomial(A, parent=P)
  end
```
In fact this now allows for two different entry points for the parent ring `P` 
of the output: First as the required `parent` keyword argument and second 
as the first argument of a method of `characteristic_polynomial` with 
an extended signature. Note that within the scope of the first method's body 
the OSCAR function `parent` is necessarily overwritten by the name of the 
keyword argument. Hence to call the actual parent of any other object, you 
must then use `Oscar.parent`. E.g. to get the `MatrixSpace` of the 
matrix `A`, write `Oscar.parent(A)`.


## Documentation

 - In general we try to follow the list of recommendations in the
   [Documentation section of the Julia manual](https://docs.julialang.org/en/v1/manual/documentation/).

 - Via the MathJax integration it is possible to use LaTeX code, and this is the
   preferred way to denote the mathematical symbols in the docstrings.


## Printing in Oscar

### The 2 + 1 print modes of Oscar
Oscar has two user print modes `detailed` and `one line` and one internal
print mode `:supercompact`. The latter is for use during recursion,
e.g. to print the `base_ring(X)` when in `one line` mode.
It exists to make sure that `one line` stays compact and human readable.

Top-level REPL printing of an object will use `detailed` mode by default
```julia
julia> X
detailed
```
Inside nested structures, e.g. inside a `Vector`, the `one line` mode is used.
```julia
julia> [X,X]
3-element Vector{TypeofX{T}}
one line
one line
one line
```

##### An Example for the 2 + 1 print modes
```
# detailed
General linear group of degree 24
  over Finite field of degree 7 over GF(29)

# one line
General linear group of degree 24 over GF(29^7)

# supercompact
General linear group
```

The print modes are specified as follows
#### Detailed printing
- the output must make sense as a standalone without context to non-specialists
- the number of output lines should fit in the terminal
- if the object is simple enough use only one line
- use indentation and (usually) `one line` to print substructures
#### One line printing
- the output must print in one line
- should make sense as a standalone without context
- variable names/generators/relations should not be printed only their number.
- Only the first word is capitalized e.g. `Polynomial ring`
- one should use `:supercompact` for nested printing in compact
- nested calls to `one line` (if you think them really necessary) should be at the end,
  so that one can read sequentially. Calls to `:supercompact` can be anywhere.
- commas must be enclosed in brackets so that printing tuples stays unambiguous
#### Super compact printing
- a user readable version of the main (mathematical) type.
- a single term or a symbol/letter mimicking mathematical notation
- should usually only depend on the type and not of the type parameters or of
  the concrete instance - exceptions of this rule are possible e.g. for `GF(2)`
- no nested printing. In particular variable names and `base_ring` must not be displayed.
  This ensures that `one line` and `:supercompact` stay compact even for complicated things.
  If you want nested printing use `one line` or `detailed`.

For further information and examples we refer you to our section [Details on
printing in Oscar](@ref).


## Deprecating functions

Sometimes it is necessary to rename a function or otherwise change it. To allow
for backwards compatibility, please then introduce a new line in the file
`src/deprecations.jl`. The syntax is as follows:
```
# Deprecated after CURRENT_RELEASE_VERSION
@deprecate old_function(args) new_function(args)
```
It is possible to transform the `args` too, if the syntax has changed. If this
process needs an auxiliary function, which otherwise is unnecessary, please add
it above:
```
# Deprecated after CURRENT_RELEASE_VERSION
function transform_args_for_new_function(args)
    # Do something
    return new_args
end
@deprecate old_function(args) new_function(transform_args_for_new_function(args))
```
The comment about the version number is only necessary if you are the first one
adding to `deprecations.jl` after a release, otherwise please add to the
existing block.

!!! note
    Please make sure to change to the new function everywhere in the existing
    OSCAR code base. Even if you think, you were the only one using the
    function, run a quick `grep` to make sure. When you are done,
    `deprecations.jl` should be the only place mentioning `old_function`. To
    make sure, you can start Julia with `--depwarn=yes` or even
    `--depwarn=error` and then run the tests.

## Approved abbreviations

- Types for rings/groups/ideals/modules/... end with `Ring`/`Group`/`Ideal`/`Module`/...
- Types for elements should have the same name as the type of the parent with
  `Elem` added;
   - Exception: `MatrixSpace` elements end with `Matrix`.
- We abbreviate certain parts of type names, according to a fixed set of
  substitutions; further abbreviations should be carefully decided upon.
- Every abbreviation must be unique; e.g. `Abs` stands for `Absolute`, and so
  must not be used for e.g. `Abstract`.
- **List of approved abbreviations**
    - absolute -> `Abs`
        - abstract -> `Abstract`
    - decorated -> `Dec`
    - group -> `Group`
    - ideal -> `Ideal`
    - localized -> `Loc`
    - matrix -> `Matrix`
    - module -> `Module`
    - multivariate polynomial -> `MPoly`
    - polynomial -> `Poly`
    - quotient -> `Quo`
    - relative -> `Rel`
    - ring ->`Ring`
    - subquotient -> `Subquo`
- If a type comes in sparse and dense variants, then call the dense type `T`
  and the sparse one `SparseT`.
