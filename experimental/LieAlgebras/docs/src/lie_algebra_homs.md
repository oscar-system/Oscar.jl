```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Lie algebra homomorphisms

Homomorphisms of Lie algebras in Oscar are represented by the type
`LieAlgebraHom`.

## Constructors

Lie algebra homomorphisms $h: L_1 \to L_2$ are constructed by providing either
the images of the basis elements of $L_1$ or a $\dim L_1 \times \dim L_2$ matrix.

```@docs
hom(::LieAlgebra{C}, ::LieAlgebra{C}, ::Vector{<:LieAlgebraElem{C}}; check::Bool=true) where {C<:FieldElem}
hom(::LieAlgebra{C}, ::LieAlgebra{C}, ::MatElem{C}; check::Bool=true) where {C<:FieldElem}
id_hom(::LieAlgebra)
zero_map(::LieAlgebra{C}, ::LieAlgebra{C}) where {C<:FieldElem}
```

## Functions

The following functions are available for `LieAlgebraHom`s:

### Basic properties
For a homomorphism $h: L_1 \to L_2$, `domain(h)` and `codomain(h)` return $L_1$ and $L_2$ respectively.

```@docs
matrix(::LieAlgebraHom{<:LieAlgebra,<:LieAlgebra{C2}}) where {C2<:FieldElem}
```

### Image
```@docs
image(::LieAlgebraHom, ::LieAlgebraElem)
image(::LieAlgebraHom)
image(::LieAlgebraHom, ::LieAlgebraIdeal)
image(::LieAlgebraHom, ::LieSubalgebra)
```

### Kernel
```@docs
kernel(::LieAlgebraHom)
```

### Composition
```@docs
compose(::LieAlgebraHom{T1,T2}, ::LieAlgebraHom{T2,T3}) where {T1<:LieAlgebra,T2<:LieAlgebra,T3<:LieAlgebra}
```

### Inverses
```@docs
is_isomorphism(::LieAlgebraHom)
inv(::LieAlgebraHom)
```
