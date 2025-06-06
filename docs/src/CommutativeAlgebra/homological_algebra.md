```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Homological Algebra

Some OSCAR functions which are fundamental to homological algebra such as the `kernel` function
for module homomorphisms and basic functions for handling chain and cochain complexes are
discussed in the module section. Building on these functions, we here introduce further OSCAR functionality
supporting computations in homological algebra.


## Pruning Modules

```@docs
prune_with_map(M::ModuleFP)
```
## Finiteness and cardinality as a set

```@docs
is_finite(M::SubquoModule{T}) where {T<:Union{ZZRingElem, FieldElem}}
```

```@docs
size(M::SubquoModule{T}) where {T<:Union{ZZRingElem, FieldElem}}
```

## Presentations

```@docs
presentation(M::ModuleFP)
```

## Representation as Cokernel

```@docs
present_as_cokernel(M::SubquoModule, task::Symbol = :none)
```

## Free Resolutions

### Types

The OSCAR type for the free resolutions discussed in this section is of parametrized form `FreeResolution{T}`.

### Constructors

```@docs
free_resolution(M::SubquoModule{T};
    length::Int = 0,
    algorithm::Symbol = T <: MPolyRingElem ? :fres : :sres) where {T <: Union{MPolyRingElem, MPolyQuoRingElem}}
```

```@docs
free_resolution(I::MPolyIdeal; length::Int = 0, algorithm::Symbol = :fres)
```

```@docs
free_resolution(Q::MPolyQuoRing; length::Int = 0, algorithm::Symbol = :fres)
```
### Data Associated to Free Resolutions

```@docs
augmented_complex(F::FreeResolution)
```

```@docs
length(F::FreeResolution)
```


## Betti Tables

Given a $\mathbb Z$-graded multivariate polynomial ring $S$, and given
a graded free resolution  with finitely generated graded free $S$-modules 

$F_i=\bigoplus_j S(-j) ^{\beta_{ij}},$

the numbers $\beta_{ij}$ are commonly known as the *graded Betti numbers*
of the resolution. A convenient way of visualizing these numbers is to write a
*Betti table* as in the example below:


```@julia
       0  1  2  3  
------------------
0    : 1  -  -  -  
1    : -  2  1  -  
2    : -  2  3  1  
------------------
total: 1  4  4  1
```

A number $i$ in the top row of the table refers to the $i$-th free 
module $F_i$ of the resolution. More precisely, the column with first 
entry $i$ lists the number of free generators
of $F_i$ in different degrees and, in the bottom row,
the total number of free generators (that is, the rank of
$F_i$). If $k$ is the first entry of a row containing 
a number $\beta$ in the column corresponding to $F_i$, 
then $F_i$ has $\beta$ generators in degree $k+i$. That is,
for a free module $F_i$ written as a direct sum as above,
$\beta$ is the number $\beta_{ij}$
with $j=k+i$. The explicit example table above indicates, for instance, 
that $F_2$ has one generator in degree 3 and three generators 
in degree 4. In total, the diagram corresponds to a 
graded free resolution of type 

$S \leftarrow S(-2)^2\oplus S(-3)^2 \leftarrow S(-3)\oplus S(-4)^3 \leftarrow S(-5) \leftarrow 0.$


```@docs
betti_table(F::FreeResolution)
```

```@docs
minimal_betti_table(F::FreeResolution{T}) where {T<:ModuleFP}
```

```@docs
minimal_betti_table(M::ModuleFP{T}) where {T<:MPolyDecRingElem}
```

## Castelnuovo-Mumford Regularity

```@docs
cm_regularity(M::ModuleFP)
```

## Homology

```@docs
homology(C::ComplexOfMorphisms{<:ModuleFP})
```

```@docs
homology(C::ComplexOfMorphisms{<:ModuleFP}, i::Int)
```

## Hom and Ext

```@docs
hom(M::ModuleFP, N::ModuleFP; algorithm::Symbol=:maps)
```

```@docs
element_to_homomorphism(f::ModuleFPElem)
```

```@docs
homomorphism_to_element(H::ModuleFP, phi::ModuleFPHom)
```

```@docs
ext(M::ModuleFP, N::ModuleFP, i::Int)
```

## Tensorproduct and Tor

```@docs
tensor_product(G::ModuleFP...; task::Symbol = :none)
```

```@docs
tor(M::ModuleFP, N::ModuleFP, i::Int)
```

## Fitting Ideals

```@docs
fitting_ideal(M::ModuleFP{T}, i::Int) where T <: MPolyRingElem
```

## Flatness

Checking flatness in OSCAR relies on characterizing flatness in terms of Fitting ideals.

```@docs
is_flat(M::ModuleFP{T}) where T <: MPolyRingElem
```

```@docs
non_flat_locus(M::ModuleFP{T}) where T <: MPolyRingElem
```

## Regular Sequence Test

```@docs
is_regular_sequence(V::Vector{T}, M::ModuleFP{T}) where T <: MPolyRingElem
```

## Koszul Complex

```@docs
koszul_matrix(V::Vector{T}, i::Int) where T <: MPolyRingElem
```

```@docs
koszul_complex(V::Vector{T}) where T <: MPolyRingElem
```

## Koszul Homology

```@docs
koszul_homology(V::Vector{T}, M::ModuleFP{T}, p::Int) where T <: MPolyRingElem
```

## Depth

The computation of depth in OSCAR relies on expressing depth in terms of  Koszul cohomology. 

```@docs
depth(I::MPolyIdeal{T}, M::ModuleFP{T}) where T <: MPolyRingElem
```






