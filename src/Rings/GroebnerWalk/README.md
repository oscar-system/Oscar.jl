# GroebnerWalk

### Usage

The following function calls perform different versions of the Gröbner Walk by computing the Gröbner basis w.r.t. the monomial order ```:degrevlex``` and converting it to the Groebner basis w.r.t. the monomial order ```:lex```.

```julia
using Singular
using Nemo
include("Examples.jl")
include("GroebnerWalk.jl")

id = katsura4()
groebnerwalk(id, :degrevlex, :lex ,:standard)
groebnerwalk(id, :degrevlex, :lex ,:pertubed, 2)
groebnerwalk(id, :degrevlex, :lex ,:fractal_combined)
groebnerwalk(id, :degrevlex, :lex ,:generic)
```

The following function call performs the standard version of the Gröbner Walk by converting the given Groebner basis w.r.t. the monomial order ``` :degrevlex``` to the Groebner basis w.r.t. the monomial order ```:lex```.

```julia
using Singular
using Nemo
include("Examples.jl")
include("GroebnerWalk.jl")

id = katsura4()
R = base_ring(id)
I = Singular.std(id, complete_reduction = true)

groebnerwalk(I, ordering_as_matrix(:degrevlex, nvars(R)), ordering_as_matrix(:lex, nvars(R)), :standard)

```

### Tests

The tests can be started with **runtests.jl**.

### Tracing the Walk

For informations about the conversions-steps use the parameter ```infoLevel```. 
