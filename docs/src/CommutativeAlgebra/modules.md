```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["modules.md"]
```

# Modules Over Polynomial Rings

In this section, the term module will always refer to a finitely presented module over a multivariate polynomial ring.
In OSCAR, the most general way of implementing such a module is that of a *subquotient*, that is,
as a submodule of a quotient of a free module. Explicitly, a *subquotient* $M$ over the ring $R$ is a module of type

$M =(\text{im } a + \text{im } b)/\text{im } b,$

where

$a:R^m ⟶R^p \;\text{ and }\; b:R^n ⟶R^p$

are two homomorphisms of free $R$-modules with the same codomain. We then refer to
- the codomain $R^p$ as the *ambient free module* of $M$,
- the images of the canonical basis vectors of $R^m$ as the *generators* of $M$, and
- the images of the canonical basis vectors of $R^n$ as the *relations* of $M$.


## Constructors

## Data Associated to Subqotients

## Elements of Subqotients

## Gröbner Bases

## Tests on Subqotients

## Basic Operations on Subqotients

## Submodules and Quotients

## Homomorphisms of Subqotients

## Operations on Homomorphisms of Subqotients

## Subquotients Related to Homomorphisms

## Presentations

## Syzygies and Free Resolutions

## Hom and Ext

## Tensorproduct and Tor

