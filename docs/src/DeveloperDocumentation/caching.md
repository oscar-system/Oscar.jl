# Caching Parent Objects

Many functions in OSCAR that construct parent objects (such as rings, modules,
groups, etc.) have an optional keyword argument `cached::Bool`. If set to
`true` then the object is put into a cache, and when the construction function
is later called again with identical inputs, then the cached object is
returned instead of creating a new object. In contrast when `cached` is set to
`false` then each time a new object is returned.

Example:
```jldoctest
julia> R1, = polynomial_ring(QQ, :x; cached = true);

julia> R2, = polynomial_ring(QQ, :x; cached = true);

julia> R1 === R2  # identical as both were created with `cached = true`
true

julia> R3, = polynomial_ring(QQ, :x; cached = false);

julia> R1 === R3  # not identical as R3 was created with `cached = false`
false

julia> R4, = polynomial_ring(QQ, :y; cached = true);

julia> R1 === R4  # not identical despite `cached = true` due to differing variable names
false
```

## Why cache parent objects?

The main reason for supporting caching of parent objects is **user convenience**:
experience shows that most mathematicians (espescially those who are not also
programmers; but it really affects all) are surprised if, say, `QQ[:x] == Q[:x]`
produces `false`.

For interactive use, it is often simply convenient: e.g. in the following example,
we use `map_coefficients` to map polynomials over the integers to polynomials
over a finite field, and the results can be added -- this is only possible because
the new polynomials have the same parent, thanks to caching.
```jldoctest
julia> Zx, x = ZZ[:x]
(Univariate polynomial ring in x over ZZ, x)

julia> F = GF(2);

julia> map_coefficients(F, x^2) + map_coefficients(F, x)
x^2 + x
```

Caching parents also has downsides. E.g. all those cached objects take up memory which
in some cases can add up to significant amounts.


## Rules for implementations

In the following we describe some rules related to caching for people implementing
parent constructor functions

1. Don't use caching in code inside OSCAR (caching is for end users!)
   - i.e., code inside OSCAR by default should always construct rings with `cached = false`.
   - In other words: internal code should not rely on caching being active.
     Usually the need for using cached parents can be overcome by allowing callers to
     pass in a parent object as an additional function argument. One may still provide a
     default value for that as a user convenience, but these default parents then should
     be created with `cached=false`.
   - Rationale: this avoids clogging the system with cached objects the user never asked
     for. It also eliminates sources of bugs: a cached ring may have attributes assigned
     that modify its behavior in a way that it is completely unexpected in code dealing
     with "newly created" ring
2. All end-user facing constructors should have a `cached::Bool` keyword argument
   with a default value, regardless of whether caching is actually supported or not.
   - if caching is supported, then `cached` should default to *true*
   - if caching is not supported, then `cached` should default to *false*
   - Rationale: this allows us to comply pro-actively with the first rule: when creating
     a parent object, you always pass in `cached = false`. If not all constructors
     support this, we can't comply with it. Even if a constructor does not support
     caching right now: this might change in the future. So by allowing the `cached`
     argument in all cases, we can write future-proof code.
3. Caches must not overflow
   - the simplest solution to achieve this is to use an `AbstractAlgebra.CacheDictType`
     instances (which really is an alias for `WeakValueDict`) together with `get_cached!`
     which automatically removes objects from caches if nothing outside the cache references
     it anymore
   - Alternatively one may offer a manual way for users to "flush" caches, but beware
     the problems this can cause when code relies on parents being cached -- yet another
     reason for rule 1.

For convenience, `Hecke` also defines these "standard rings" for use in functions
like `cyclotomic_polynomial`
```
module Globals
  using Hecke
  const Qx, _ = polynomial_ring(QQ, :x, cached = false)
  const Zx, _ = polynomial_ring(ZZ, :x, cached = false)
  const Zxy, _ = polynomial_ring(ZZ, [:x, :y], cached = false)
end
```
You can use these in your own code as well, or imitate this pattern if convenient.

As always, if in doubt what to do, please ask.
