```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Enumeration of isometries

One of the main features of this project is the enumeration of even lattices
equipped with an isometry of finite order. It is based on an implementation of
the algorithms of [BH23](@cite) for isometries whose order has at most 2 prime
divisors. Note that such algorithms have been generalized to cover all finite
orders using an iterative procedure and Nikulin's theory of equivariant
primitive extensions (also known as _orthogonal primitive extensions_).

We guide the user here to the global aspects of the available theory, and we
refer to the paper [BH23](@cite) for further reference.

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

We use Definition 4.13. and Algorithm 1 of [BH23](@cite) to implement the
necessary tools for working with admissible triples. Most of the computations
consists of local genus symbol manipulations and combinatorics. The code also
relies on enumeration of integer genera with given signatures, determinant
and bounded scale valuations for the Jordan components at all the relevant
primes (see [`integer_genera`](@ref)).

```@docs
admissible_triples(::ZZGenus, ::Int)
is_admissible_triple(::ZZGenus, ::ZZGenus, ::ZZGenus, ::Int)
```

Note that admissible triples are mainly used for enumerating lattices with
isometry of a given order and in a given genus.

## Enumeration functions

We give an overview of the functions implemented for the enumeration of the
isometries of even $\mathbb{Z}$-lattices. For more details such as the proof
of the algorithms and the theory behind them, we refer to [BH23](@cite) and
references therein.

### The hermitian case

For an irreducible reciprocal polynomial $\chi$ and a genus symbol $G$
of integral integer lattices, if the equation order $\mathbb{Z}[\chi]$ is maximal,
one can compute representatives of isomorphism classes of lattices with isometry
$(L, f)$ such that $L\in G$ and $\chi(f) = 0$.

```@docs
representatives_of_hermitian_type(::Union{ZZLat, ZZGenus}, ::Union{ZZPolyRingElem, QQPolyRingElem}, ::Int)
```

In the case of finite order isometries, when $\chi$ is cyclotomic, one can use
as a shortcut the following function instead:

```@docs
representatives_of_hermitian_type(::Union{ZZGenus, ZZLat}, ::Int, ::Int)
```

### The generic case

The algorithms from [BH23](@cite) are specialized on the requirement for the
input and regular users are not expected to known which function to choose for
their purpose. We therefore provide a generic function which allows one to
compute a complete of representatives for the isomorphism classes of even
lattices equipped with an isometry of finite order.

```@docs
enumerate_classes_of_lattices_with_isometry(::Union{ZZGenus, ZZLat}, ::Int)
```

As a remark: if $n = p_1^{e_1}p_2^{e_2}...p_k^{e_k}$ is the chosen order, with
$p_1 < p_2 < ... < p_k$ distinct prime numbers and $e_i > 0$ for all $1\leq
i\leq k$, then the previous function computes first iteratively representatives
for all classes in the given genus with an isometry of order finite
$p_k^{e_k}$. Then, the function iteratively increases the order to be
$p_{k-1}^{e_{k-1}}p_k^{e_k}$, and so on until $n$.

In particular, the current naive approach requires to determine representatives
for the isomorphism classes of even lattices with isometry of order dividing
$n$, for some divisors of $n$. Hence, if $n$ has many divisors, such
computations can be expensive. We do not recommend to use such a method for
lattices of large rank or isometries of high order as the computations could
possibly take several weeks. Interested users may look into the next section
for more advanced methods, to be used on a larger framework.

## Underlying machinery

Here is a list of the algorithmic machinery provided by Brandhorst and Hofmann
in [BH23](@cite), and further extended by Muller, which is used by the previous
enumerating functions. The functions goes from the more specialized method, to
the most generic one.

```@docs
representatives_of_hermitian_type(::ZZLatWithIsom, ::Int, ::Int)
splitting_of_hermitian_type(::ZZLatWithIsom, ::Int, ::Int)
splitting_of_prime_power(::ZZLatWithIsom, ::Int, ::Int)
splitting_of_pure_mixed_prime_power(::ZZLatWithIsom, ::Int)
splitting_of_mixed_prime_power(::ZZLatWithIsom, ::Int, ::Int)
splitting(::ZZLatWithIsom, ::Int, ::Int)
```

Note that an important feature from the theory in [BH23](@cite) is the notion
of *admissible gluings* and equivariant primitive embeddings for admissible
triples. In the next chapter, we present the methods regarding Nikulins's
theory on primitive embeddings and their equivariant version. We use this
basis to introduce the method
[`admissible_equivariant_primitive_extensions`](@ref) (Algorithm 2 in
[BH23](@cite)) which is the major tool making the previous enumeration
possible and fast, from an algorithmic point of view.

