```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

# Lie algebra module homomorphisms

Homomorphisms of Lie algebra modules in Oscar are represented by the type
`LieAlgebraModuleHom`.

## Constructors

Homomorphisms of modules over the same Lie algebra
$h: V_1 \to V_2$ are constructed by providing either the images of the basis elements of $V_1$
or a $\dim V_1 \times \dim V_2$ matrix.

```@docs
hom(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::Vector{<:LieAlgebraModuleElem{C}}; check::Bool=true) where {C<:RingElement}
hom(::LieAlgebraModule{C}, ::LieAlgebraModule{C}, ::MatElem{C}; check::Bool=true) where {C<:RingElement}
identity_map(::LieAlgebraModule)
```

## Functions

The following functions are available for `LieAlgebraModuleHom`s:

### Basic properties
For a homomorphism $h: V_1 \to V_2$, `domain(h)` and `codomain(h)` return $V_1$ and $V_2$ respectively.

```@docs
matrix(::LieAlgebraModuleHom{<:LieAlgebraModule,<:LieAlgebraModule{C2}}) where {C2<:RingElement}
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
