```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_cones.md"]
```

# Cones

Test

## Introduction

Let $\mathbb{F}$ be an ordered field; the most prominent case here is where $\mathbb{F}=\mathbb{Q}$ are the rational numbers.

A set $C \subseteq \mathbb{F}^n$ is called a (polyhedral) cone if it can be written as the set of nonnegative linear combinations of finitely many vectors in $\mathbb{F}^n$.
Eqivalently, cones can be written as the intersection of finitely many homogeneous linear inequalities.

Any cone is a special case of a polyhedron.
Conversely, intersecting a cone with a suitable affine hyperplane yields a polyhedron whose faces are in bijection swith the faces of the cone.
Going back and forth between polyhedra and their homogenizations, the cones, is a frequent operation.
This is one reason for keeping cones as a distinct type.

## Construction
