```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Tropical polyhedra

## Introduction
A tropical polytope $\mathrm{tconv}(V)$ is the tropical convex hull of finitely many points $V$ which is defined as the set of
tropical linear combinations of the points $v_i\in V$ with coefficients in $\mathbb{R}$.
The point set may be taken from either the tropical projective torus $\mathbb{R}^d/\mathbb{R}\mathbf{1}$ 
or the tropical projective space $$\mathbb{TP}^{d-1} := (\mathbb{T}^d\setminus\{\mathbf{\infty}\})/\mathbb{R}\mathbf{1}.$$
Likewise, $\mathrm{tconv}(V)$ may be considered in either $\mathbb{R}^d/\mathbb{R}\mathbf{1}$ or $\mathbb{TP}^{d-1}$.

For more details see
- Chapter 5.2 in [MS15](@cite)
- Chapter 5 in [Jos21](@cite)

#### Note:
- We offer two types with slightly different semantics, `TropicalPolyhedron` and `TropicalPointConfiguration`.
  - The `TropicalPointConfiguration` type represents the information of the covector decomposition of the tropical projective torus
    induced by $V$.
  - A `TropicalPolyhedron` refers to the actual tropical convex hull of the point set $V$ which carries an induced covector decomposition,
    which is a polyhedral subcomplex of the covector decomposition for the corresponding `TropicalPointConfiguration`.
- Currently, OSCAR supports
  - point sets $V$ with entirely finite coordinates
  - point sets $V$ with infinite coordinates and $\mathrm{tconv}(V)$ inside the tropical projective torus
- Functionality for tropical polytopes over the entire tropical projective space is planned for a future version

## Constructors
Both `TropicalPolyhedron` and `TropicalPointConfiguration` can be constructed from the following data:
- a rational matrix and specifying min- or max-convention
- a matrix over a tropical semiring
- a vector of points with coordinates in a tropical semiring

Additionally, one can convert between `TropicalPolyhedron` and `TropicalPointConfiguration`. This way, 
the underlying low-level object in `polymake` and already computed information is preserved.
```@docs
tropical_convex_hull
tropical_point_configuration
```

## Properties
`TropicalPolyhedron` and `TropicalPointConfiguration` are closely related as being given by a finite set of points, 
but refer to different interpretations of this point configuration. In general, the properties of a `TropicalPolyhedron`
refer to the tropical convex hull of these points while `TropicalPointConfiguration` refers to the 
covector decomposition of the entire tropical projective torus.

```@docs
ambient_dim(P::Union{TropicalPolyhedron,TropicalPointConfiguration})
```

### Tropical polyhedra
```@docs
vertices(as::Type{PointVector{T}}, P::TropicalPolyhedron) where {T<:TropicalSemiringElem}
n_vertices(P::TropicalPolyhedron)
dim(P::TropicalPolyhedron)
maximal_covectors(as::Type{IncidenceMatrix}, P::TropicalPolyhedron)
covector_decomposition(P::TropicalPolyhedron)
is_bounded(P::TropicalPolyhedron)
```

### Tropical point configurations
```@docs
points(P::TropicalPointConfiguration{M}) where {M<:MinOrMax}
n_points(P::TropicalPointConfiguration)
maximal_covectors(as::Type{IncidenceMatrix}, P::TropicalPointConfiguration)
covector_decomposition(P::TropicalPointConfiguration)
```
