```@meta
CurrentModule = Oscar
```

One of the main feature of this project is the enumeration of lattices with
isometry of finite order with at most two prime divisors. This is the content
of [BH23] which has been implemented. We guide the user here to the global
aspects of the available theory, and we refer to the paper [BH23] for further
reference.

# Admissible triples

Roughly speaking, for a prime number $p$, a *$p$-admissible triple* `(A, B, C)`
is a triple of integer lattices such that,in certain cases, `C` can be obtained
has a primitive extension $A \perp B \to C$ where one can glue along
$p$-elementary subgroups of the respective discriminant groups of `A` and `B`.
Note that not all admissible triples satisfy this extension property.

For instance, if $f$ is an isometry of an integer lattice `C` of prime order
`p`, then for $A := \Ker \Phi_1(f)$ and $B := \Ker \Phi_p(f)$, one has that
`(A, B, C)` is $p$-admissible (see [BH13, Lemma 4.15.]).

We say that a triple `(AA, BB, CC)` of genus symbols for integer lattices is
*$p$-admissible* if there are some lattices $A \in AA$, $B \in BB$ and
$C \in CC$ such that $(A, B, C)$ is $p$-admissible.

We use Definition 4.13. and Algorithm 1 of [BH23] to implement the necessary
tools for working with admissible triples. Most of the computations consists of
local genus symbol manipulations and combinatorics. The code also relies on
enumeration of integer genera with given signatures, determinant and bounded
scale valuations for the Jordan components at all the relevant primes (see
[`integer_genera`](@ref)).

```@docs
admissible_triples(:::ZZGenus, p::Integer)
is_admissible_triple(::ZZGenus, ::ZZGenus, ::ZZGenus, ::Integer)
```

Note that admissible triples are mainly used for enumerating lattices with
isometry of a given order and in a given genus.

# Enumeration functions

We give an overview of the functions implemented for the enumeration of the
isometries of integral integer lattices. For more details such as the proof of
the algorithms and the theory behind them, we refer to the reference paper [BH23].

## Global function

As we will see later, the algorithms from [BH23] are specialized on the
requirement for the input and regular users might not be able to choose between
the functions available. We therefore provide a general function which
allows one to enumerate lattices with isometry of a given order and in a given
genus. The only requirements are to provide a genus symbol, or a lattice from
this genus, and the order wanted (as long as the number of distinct prime
divisors is at most 2).

```@docs
enumerate_classes_of_lattices_with_isometry(::ZZLat, ::Hecke.IntegerUnion)
```

As a remark: if $n = p^dq^e$ is the chosen order, with $p < q$ prime numbers,
the previous function computes first iteratively representatives for all classes
with isometry in the given genus of order $q^e$. Then, the function increases
iteratively the order up to $p^dq^e$.

## Underlying machinery

Here is a list of the algorithmic machinery provided by [BH23] used previously
to enumerate lattices with isometry:

```@docs
representatives_of_hermitian_type(::ZZLatWithIsom, ::Int)
splitting_of_hermitian_prime_power(::ZZLatWithIsom, ::Int)
splitting_of_prime_power(::ZZLatWithIsom, ::Int, ::Int)
splitting_of_pure_mixed_prime_power(::ZZLatWithIsom, ::Int)
splitting_of_mixed_prime_power(::ZZLatWithIsom, ::Int, ::Int)
```

Note that an important feature from the theory in [BH23] is the notion of
*admissible gluings* and equivariant primitive embeddings for admissible triples.
In the next chapter, we will develop about the features regarding primitive
embeddings and their equivariant version. We use this basis to introduce the
method [`admissible_equivariant_primitive_extension`](@ref) (Algorithm 2 in
[BH23]) which is the major tool making the previous enumeration possible and
fast, from an algorithmic point of view.
