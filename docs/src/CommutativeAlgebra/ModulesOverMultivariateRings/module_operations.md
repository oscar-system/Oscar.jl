```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["module_operations.md"]
```

# Operations on Modules

## Subquotients Related to Homomorphisms

### Kernel

```@docs
kernel(a::ModuleFPHom)
```

### Image

```@docs
image(a::ModuleFPHom)
```

### Cokernel

```@docs
cokernel(a::ModuleFPHom)
```

## Chain and Cochain Complexes

The general OSCAR type `ChainComplex{T}` allows one to model both chain complexes and cochain complexes
(the `T` refers to the type of the differentials of the complex). In the context of commutative algebra, we handle
complexes of modules and module homomorphisms over multivariate polynomial rings.
In this section, we first show how to create such complexes. Then we discuss functionality for dealing with the constructed
complexes, mainly focusing on chain complexes. Cochain complexes can be handled similarly.

### Constructors

```@docs
chain_complex(V::ModuleFPHom...; start::Int = 0)
```

```@docs
cochain_complex(V::ModuleFPHom...; start::Int = 0)
```

### Data Associated to Chain Complexes

Given a chain complex `C`,
- `range(C)` refers to the range of `C`,
- `obj(C, i)` and `C[i]` to the `i`-th module of `C`, and
- `map(C, i)` to the `i`-th differential of `C`.

#### Examples

```@repl oscar
R, (x,) = PolynomialRing(QQ, ["x"]);
F = free_module(R, 1);
A, _ = quo(F, [x^4*F[1]]);
B, _ = quo(F, [x^3*F[1]]);
a = hom(A, B, [x^2*B[1]]);
b = hom(B, B, [x^2*B[1]]);
C = chain_complex([a, b]; start =3);
range(C)
C[5]
map(C, 5)
```

### Operations on Chain Complexes

```@julia
shift(C::ChainComplex{T}, n::Int) where T
```

Return the complex obtained from `C` by shifting the homological degrees `d` steps.

#### Examples

```@repl oscar
R, (x,) = PolynomialRing(QQ, ["x"]);
F = free_module(R, 1);
A, _ = quo(F, [x^4*F[1]]);
B, _ = quo(F, [x^3*F[1]]);
a = hom(A, B, [x^2*B[1]]);
b = hom(B, B, [x^2*B[1]]);
C = chain_complex([a, b]; start = 3);
range(C)
D = shift(C, 3);
range(D)
```

```@docs
hom(C::ChainComplex{ModuleFP}, M::ModuleFP)
```

```@docs
hom(M::ModuleFP, C::ChainComplex{ModuleFP})
```

```@docs
tensor_product(C::ChainComplex{ModuleFP}, M::ModuleFP)
```

```@docs
tensor_product(M::ModuleFP, C::ChainComplex{ModuleFP})
```

### Tests on Chain Complexes

### Maps of Chain Complexes

#### Types

#### Constructors

## Direct Sums and Products

```@docs
direct_sum(M::ModuleFP{T}...; task::Symbol = :sum) where T
```

```@docs
direct_product(M::ModuleFP{T}...; task::Symbol = :prod) where T
```

## Presentations

```@docs
presentation(M::ModuleFP)
```

## Syzygies and Free Resolutions

```@docs
free_resolution(M::SubQuo; ordering::ModuleOrdering = default_ordering(M),
    length::Int=0, algorithm::Symbol=:fres)
```

## Homology

```@docs
homology(C::ChainComplex{<:ModuleFP})
```

```@docs
homology(C::ChainComplex{<:ModuleFP}, i::Int)
```

## Hom and Ext

```@docs
hom(M::ModuleFP, N::ModuleFP, alg::Symbol=:maps)
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

