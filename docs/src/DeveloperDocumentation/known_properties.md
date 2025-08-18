```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```
# Known properties

Many high-level objects in Oscar such as ideals in polynomial rings, algebraic
schemes, etc. have properties which can be computed and possibly be cached.
For example, `krull_dim(I::MPolyIdeal)` returns the Krull dimension of `I`
and, if not already present, caches the result in order to avoid a second
computation of this property in the future.

In various settings it is then convenient to assess whether or not a certain
property is already known, or whether calling it would trigger a non-trivial
computation. As an example, think of a user-facing method to compute the
normalization of the zero locus of an ideal `I` as above. There are
specialized routines for curves which are faster than the generic ones for
affine schemes of arbitrary dimension. Thus if the dimension of `I` was
already known to be one, one could deflect for such a special routine
directly. Yet, computation of an ideal's Krull dimension involves a Gröbner
basis computation which, if not already done, may clog the main process for an
indefinite amount of time. This risk of blockage only for reasons of algorithm
selection might be considered too high: The user "only" wants a
normalization, but Oscar gets stuck on gathering information on which
implementation to choose from. For the user's convenience, one could then
write
```julia
  function my_normalization(I::MPolyIdeal)
    if is_known(krull_dim, I) && is_one(krull_dim(I))
      return special_and_fast_normalization_for_curves(I)
    else
      return generic_but_slow_normalization(I)
    end
  end
```
using methods of `is_known` for this setup.
```@docs
    is_known
```
Some further remarks:

* The `is_known` methods have to be implemented manually, and do not
  automatically exist for all arguments.
* The above example does not adhere to the style guide in Oscar for reasons
  different than the purpose illustrated here. For instance, the output does not
  depend on the mathematical content of the input alone, but on the status of
  cached objects. This can lead to hard-to-debug code down the road and should,
  in general, be avoided. 
* Other systems which are already integrated in Oscar have their own caching
  mechanisms; for instance the attributes system in GAP, or the storage of
  attributes in Polymake objects. The function `is_known` is a priori not tied
  to any of these! However, if needed, special methods of `is_known` can be
  implemented which use such internals. 
* The function `is_known` is not user-facing. It provides an Oscar internal
  basis to assess known properties so that programmers can extend this
  functionality for their purposes if needed, using one and the same formulation
  throughout the Oscar system. In particular, we do not require `is_known` to be
  supported for any property of any object and the user can not expect this. 
  
