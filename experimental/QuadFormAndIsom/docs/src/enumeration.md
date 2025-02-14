```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Enumeration of isometries

One of the main features of this project is the enumeration of even lattices with
isometry of finite order with at most two prime divisors. This is the content
of [BH23](@cite) which has been implemented. We guide the user here to the global
aspects of the available theory, and we refer to the paper [BH23](@cite) for further
reference.

## Admissible triples

Roughly speaking, for a prime number $p$, a *$p$-admissible triple* $(A, B, C)$
is a triple of integer lattices such that, in some cases, $C$ can be obtained
as a primitive extension $A \oplus B \to C$ where one can glue along
$p$-elementary subgroups of the respective discriminant groups of $A$ and $B$.
Note that not all admissible triples satisfy this extension property.

For instance, if $f$ is an isometry of an integer lattice $C$ of prime order
$p$, then for $A := \ker \Phi_1(f)$ and $B := \ker \Phi_p(f)$, one has that
$(A, B, C)$ is $p$-admissible (see Lemma 4.15. in [BH23](@cite)).

We say that a triple $(G_A, G_B, G_C)$ of genus symbols for integer lattices is
*$p$-admissible* if there are some lattices $A \in G_A$, $B \in G_B$ and
$C \in G_C$ such that $(A, B, C)$ is $p$-admissible.

We use Definition 4.13. and Algorithm 1 of [BH23](@cite) to implement the necessary
tools for working with admissible triples. Most of the computations consists of
local genus symbol manipulations and combinatorics. The code also relies on
enumeration of integer genera with given signatures, determinant and bounded
scale valuations for the Jordan components at all the relevant primes (see
[`integer_genera`](@ref)).

```@docs
admissible_triples(::ZZGenus, p::Integer)
is_admissible_triple(::ZZGenus, ::ZZGenus, ::ZZGenus, ::Integer)
```

Note that admissible triples are mainly used for enumerating lattices with
isometry of a given order and in a given genus.

## Enumeration functions

We give an overview of the functions implemented for the enumeration of the
isometries of integral integer lattices. For more details such as the proof of
the algorithms and the theory behind them, we refer to the reference paper
[BH23](@cite).

### The hermitian case

For an irreducible reciprocal polynomial $\chi$ and a genus symbol $G$
of integral integer lattices, if the equation order $\mathbb{Z}[\chi]$ is maximal,
one can compute representatives of isomorphism classes of lattices with isometry
$(L, f)$ such that $L\in G$ and $\chi(f) = 0$.

```@docs
representatives_of_hermitian_type(::Union{ZZLat, ZZGenus}, ::Union{ZZPolyRingElem, QQPolyRingElem}, ::Bool)
```

In the case of finite order isometries, when $\chi$ is cyclotomic, one can use
as a shortcut the following function instead:

```@docs
representatives_of_hermitian_type(::Union{ZZGenus, ZZLat}, ::Int, ::Bool)
```

### Orders with two prime divisors

As we will see later, the algorithms from [BH23](@cite) are specialized on the
requirement for the input and regular users might not be able to choose between
the functions available. We therefore provide a general function which
allows one to enumerate lattices with isometry of a given order and in a given
genus. The only requirements are to provide a genus symbol, or a lattice from
this genus, and the order wanted (as long as the number of distinct prime
divisors is at most 2).

```@docs
enumerate_classes_of_lattices_with_isometry(::ZZLat, ::IntegerUnion)
enumerate_prime_power_isometries(::ZZLat, ::IntegerUnion, ::IntegerUnion)
```

As a remark: if $n = p^dq^e$ is the chosen order, with $p < q$ prime numbers,
the previous function computes first iteratively representatives for all classes
with isometry in the given genus of order $q^e$. Then, the function increases
iteratively the order up to $p^dq^e$.

### Underlying machinery

Here is a list of the algorithmic machinery provided by [BH23](@cite) used
previously to enumerate lattices with isometry. The following correspond
respectively to Algorithms 3, 4, 5, 6 and 7 of the aforementioned paper:

```@docs
representatives_of_hermitian_type(::ZZLatWithIsom, ::Int)
splitting_of_hermitian_prime_power(::ZZLatWithIsom, ::Int)
splitting_of_prime_power(::ZZLatWithIsom, ::Int, ::Int)
splitting_of_pure_mixed_prime_power(::ZZLatWithIsom, ::Int)
splitting_of_mixed_prime_power(::ZZLatWithIsom, ::Int, ::Int)
```

Note that an important feature from the theory in [BH23](@cite) is the notion of
*admissible gluings* and equivariant primitive embeddings for admissible triples.
In the next chapter, we present the methods regarding Nikulins's theory on primitive
embeddings and their equivariant version. We use this basis to introduce the
method [`admissible_equivariant_primitive_extensions`](@ref) (Algorithm 2 in
[BH23](@cite)) which is the major tool making the previous enumeration
possible and fast, from an algorithmic point of view.

