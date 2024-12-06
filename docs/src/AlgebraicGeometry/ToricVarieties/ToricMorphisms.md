```@meta
CurrentModule = Oscar
```

# ToricMorphisms

A class of morphisms among toric varieties are described by certain lattice morphisms.
Let $N_1$ and $N_2$ be lattices and $\Sigma_1$, $\Sigma_2$ fans in $N_1$ and $N_2$
respectively. A $\mathbb{Z}$-linear map

$$\overline{\phi} \colon N_1 \to N_2$$

is said to be compatible with the fans $\Sigma_1$ and $\Sigma_2$ if for every cone
$\sigma_1 \in \Sigma_1$, there exists a cone $\sigma_2 \in \Sigma_2$ such that
$\overline{\phi}_{\mathbb{R}}(\sigma_1) \subseteq \sigma_2$.

By theorem 3.3.4 [CLS11](@cite), such a map $\overline{\phi}$ induces a morphism
$\phi \colon X_{\Sigma_1} \to X_{\Sigma_2}$ of the toric varieties, and those
morphisms are exactly the toric morphisms.


## Constructors


### Generic constructors with specified codomain

```@docs
toric_morphism(domain::NormalToricVarietyType, mapping_matrix::ZZMatrix, codomain::NormalToricVarietyType; check=true)
toric_morphism(domain::NormalToricVarietyType, grid_morphism::FinGenAbGroupHom, codomain::NormalToricVarietyType; check=true)
```

### Special constructors

```@docs
toric_identity_morphism(variety::NormalToricVarietyType)
```


## Attributes of Toric Morphisms

### General attributes

```@docs
domain(tm::ToricMorphism)
image(tm::ToricMorphism)
codomain(tm::ToricMorphism)
grid_morphism(tm::ToricMorphism)
morphism_on_torusinvariant_weil_divisor_group(tm::ToricMorphism)
morphism_on_torusinvariant_cartier_divisor_group(tm::ToricMorphism)
morphism_on_class_group(tm::ToricMorphism)
morphism_on_picard_group(tm::ToricMorphism)
covering_morphism(f::ToricMorphism)
```

### Special attributes of toric varieties

To every toric variety $v$ we can associate a special toric variety, the
Cox variety. By definition, the Cox variety is such that the mapping matrix of
the toric morphism from the Cox variety to the variety $v$ is simply
given by the ray generators of the variety $v$. Put differently,
if there are exactly $N$ ray generators for the fan of $v$, then the
Cox variety of $v$ has a fan for which the ray generators are the standard basis
of $\mathbb{R}^N$ and the maximal cones are one to one to the maximal cones of
the fan of $v$.

```@docs
morphism_from_cox_variety(variety::NormalToricVarietyType)
cox_variety(variety::NormalToricVarietyType)
```
