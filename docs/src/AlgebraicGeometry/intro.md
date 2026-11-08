# Introduction

The algebraic geometry part of OSCAR provides functionality for dealing with

* affine algebraic sets and varieties
* projective algebraic sets and varieties
* schemes
* toric varieties,
* toric schemes,

Computations in affine and projective algebraic geometry rely on [Commutative Algebra](@ref).


Similarly, most algorithms for toric varieties and schemes are are based on
[Polyhedral Geometry](@ref).

General textbooks offering details on the theory of varieties and schemes include:
- [Har77](@cite)
- [The Stacks Project](https://stacks.math.columbia.edu)


## Tutorials

We encourage you to take a look at our [tutorials](https://www.oscar-system.org/tutorials/). For instance,
there is a [tutorial for toric geometry in OSCAR](https://www.oscar-system.org/tutorials/ToricGeometry/).


## Conventions

### Projectivization

There are two opposite conventions in common use when defining the projectivization. For details, see [https://stacks.math.columbia.edu/tag/01OA](https://stacks.math.columbia.edu/tag/01OA) (search for "Warning").

For an example, look at proposition 3 in [https://arxiv.org/abs/1501.04049](https://arxiv.org/abs/1501.04049). This proposition states that any elliptically fibred K3-surface can be described as hypersurface in the space $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}(0) \oplus \mathcal{O}_{\mathbb{P}^1}(-4) \oplus \mathcal{O}_{\mathbb{P}^1}(-6))$. Authors, that apply the opposite convention, would say that any elliptically fibred K3-surface is a hypersurface in $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}(0) \oplus \mathcal{O}_{\mathbb{P}^1}(4) \oplus \mathcal{O}_{\mathbb{P}^1}(6))$. Note the opposite signs.

Irrespective of the signs used to denote the ambient space, there is always agreement on the hypersurface equation and the transformation behavior of the local coordinates. In the above example, the hypersurface equation is

$$z y^2 = x^3 + \alpha(s,t) x z^2 + \beta(s,t) z^3 \, ,$$

where $(s,t)$ are coordinates on $\mathbb{P}^1$ and $\alpha(s,t)$, $\beta(s,t)$ are homogeneous polynomials in $s$, $t$ of degrees 8 and 12, respectively. The projective coordinates $[x : y : z]$ of the ambient space transform as sections of the bundles $\mathcal{O}_{\mathbb{P}^1}(4)$, $\mathcal{O}_{\mathbb{P}^1}(6)$ and $\mathcal{O}_{\mathbb{P}^1}(0)$, respectively.
