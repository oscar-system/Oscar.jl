```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["intro.md"]
```

# Introduction

In this chapter, we introduce structures and functionality for algebraic geometry. Among others, we support:
* schemes,
* toric varieties,
* toric schemes,
* tropical geometry.

We refer to the individual sections for more details.


## Conventions

### Projectivization

There are two opposite conventions in common use when defining the projectivization. For details, see [https://stacks.math.columbia.edu/tag/01OA](https://stacks.math.columbia.edu/tag/01OA) (search for "Warning").

For an example, look at proposition 3 in [https://arxiv.org/abs/1501.04049](https://arxiv.org/abs/1501.04049). This proposition states that any elliptically fibred K3-surface can be described as hypersurface in the space $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}(0) \oplus \mathcal{O}_{\mathbb{P}^1}(-4) \oplus \mathcal{O}_{\mathbb{P}^1}(-6))$. Authors, that apply the opposite convention, would say that any elliptically fibred K3-surface is a hypersurface in $\mathbb{P}(\mathcal{O}_{\mathbb{P}^1}(0) \oplus \mathcal{O}_{\mathbb{P}^1}(4) \oplus \mathcal{O}_{\mathbb{P}^1}(6))$. Note the opposite signs.

Irrespective of the signs used to denote the ambient space, there is always agreement on the hypersurface equation and the transformation behavior of the local coordinates. In the above example, the hypersurface equation is

$$z y^2 = x^3 + \alpha(s,t) x z^2 + \beta(s,t) z^3 \, ,$$

where $(s,t)$ are coordinates on $\mathbb{P}^1$ and $\alpha(s,t)$, $\beta(s,t)$ are homogeneous polynomials in $s$, $t$ of degrees 8 and 12, respectively. The projective coordinates $[x : y : z]$ of the ambient space transform as sections of the bundles $\mathcal{O}_{\mathbb{P}^1}(4)$, $\mathcal{O}_{\mathbb{P}^1}(6)$ and $\mathcal{O}_{\mathbb{P}^1}(0)$, respectively.
