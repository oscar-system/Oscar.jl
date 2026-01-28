```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Lie algebra module homomorphisms

Homomorphisms of Lie algebra modules in Oscar are represented by the type
`LieAlgebraModuleHom`.

## Constructors

Homomorphisms of modules over the same Lie algebra
$h: V_1 \to V_2$ are constructed by providing either the images of the basis elements of $V_1$
or a $\dim V_1 \times \dim V_2$ matrix.

```@docs
hom(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::Vector{<:LieAlgebraModuleElem{C}}; check::Bool=true) where {C<:FieldElem}
hom(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::MatElem{C}; check::Bool=true) where {C<:FieldElem}
id_hom(::LieAlgebraModule)
zero_map(::LieAlgebraModule{C}, ::LieAlgebraModule{C}) where {C<:FieldElem}
```

## Functions

The following functions are available for `LieAlgebraModuleHom`s:

### Basic properties
For a homomorphism $h: V_1 \to V_2$, `domain(h)` and `codomain(h)` return $V_1$ and $V_2$ respectively.

```@docs
matrix(::LieAlgebraModuleHom{<:LieAlgebraModule,<:LieAlgebraModule{C2}}) where {C2<:FieldElem}
```

### Image
```@docs
image(::LieAlgebraModuleHom, ::LieAlgebraModuleElem)
```

### Composition
```@docs
compose(::LieAlgebraModuleHom{T1,T2}, ::LieAlgebraModuleHom{T2,T3}) where {T1<:LieAlgebraModule,T2<:LieAlgebraModule,T3<:LieAlgebraModule}
```

### Inverses
```@docs
is_isomorphism(::LieAlgebraModuleHom)
inv(::LieAlgebraModuleHom)
```

### Hom constructions
Lie algebra module homomorphisms support `+` and `-` if they have the same domain and codomain.

```@docs
canonical_injections(::LieAlgebraModule)
canonical_injection(::LieAlgebraModule, ::Int)
canonical_projections(::LieAlgebraModule)
canonical_projection(::LieAlgebraModule, ::Int)
hom_direct_sum(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::Matrix{<:LieAlgebraModuleHom}) where {C<:FieldElem}
hom_tensor(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::Vector{<:LieAlgebraModuleHom}) where {C<:FieldElem}
hom(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::LieAlgebraModuleHom) where {C<:FieldElem}
```
