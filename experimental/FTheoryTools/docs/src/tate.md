```@meta
CurrentModule = Oscar
```

# Global Tate models

## Introduction

A global Tate model describes a particular form of an elliptic fibration.
We focus on an elliptic fibration over a base ``B``. Consider
the weighted projective space ``\mathbb{P}^{2,3,1}`` with coordinates
``x, y, z``. In addition, consider
* ``a_1 \in H^0( B_3, \overline{K}_{B} )``,
* ``a_2 \in H^0( B_3, \overline{K}_{B}^{\otimes 2} )``,
* ``a_3 \in H^0( B_3, \overline{K}_{B}^{\otimes 3} )``,
* ``a_4 \in H^0( B_3, \overline{K}_{B}^{\otimes 4} )``,
* ``a_6 \in H^0( B_3, \overline{K}_{B}^{\otimes 6} )``.
Then form a ``\mathbb{P}^{2,3,1}``-bundle over ``B`` such that
* ``x`` transforms as a section of ``2 \overline{K}_{B}``,
* ``y`` transforms as a section of ``3 \overline{K}_{B}``,
* ``z`` transforms as a section of ``0 \overline{K}_{B} = \mathcal{O}_{B}``.
In this 5-fold ambient space, a global Tate model is the hypersurface defined
by the vanishing of the Tate polynomial
``P_T = x^3 - y^2 - x y z a_1 + x^2 z^2 a_2 - y z^3 a_3 + x z^4 a_4 + z^6 a_6``.

Crucially, for non-trivial F-theory settings, the elliptic fibration in question must
be singular. In fact, by construction, one usually engineers certain singularities.
For this, vanishing orders of the sections ``a_i`` above need to specified. The
following table---often referred to as the Tate table and taken from
[Wei10](@cite)---summarizes the singularities introduced by certain vanishing orders:

| sing. type | ``\mathrm{ord}(\Delta)`` | singularity | group ``G`` | ``a_1`` | ``a_2`` | ``a_3`` | ``a_4`` | ``a_6`` |
| :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- | :----------- |
| ``I_0`` | ``0`` |  |  | ``0`` | ``0`` | ``0`` | ``0`` | ``0`` |
| ``I_1`` | ``1`` |  |  | ``0`` | ``0`` | ``1`` | ``1`` | ``1`` |
| ``I_2`` | ``2`` | ``A_1`` | ``SU(2)`` | ``0`` | ``0`` | ``1`` | ``1`` | ``2`` |
| ``I_{2k}^{ns}`` | ``2k`` | ``C_k`` | ``Sp(k)`` | ``0`` | ``0`` | ``k`` | ``k`` | ``2k`` |
| ``I_{2k}^s`` | ``2k`` | ``A_{2k-1}`` | ``SU(2k)`` | ``0`` | ``1`` | ``k`` | ``k`` | ``2k`` |
| ``I_{2k+1}^{ns}`` | ``2k+1`` | | ``Sp(k)`` | ``0`` | ``0`` | ``k+1`` | ``k+1`` | ``2k+1`` |
| ``I_{2k+1}^{s}`` | ``2k+1`` | ``A_{2k}`` | ``SU(2k+1)`` | ``0`` | ``1`` | ``k`` | ``k+1`` | ``2k+1`` |
| ``II`` | ``2`` | | | ``1`` | ``1`` | ``1`` | ``1`` | ``1`` |
| ``III`` | ``3`` | ``A_1`` | ``SU(2)`` | ``1`` | ``1`` | ``1`` | ``1`` | ``2`` |
| ``IV^{ns}`` | ``4`` |  | ``Sp(1)`` | ``1`` | ``1`` | ``1`` | ``2`` | ``2`` |
| ``IV^s`` | ``4`` | ``A_2`` | ``SU(3)`` | ``1`` | ``1`` | ``1`` | ``2`` | ``3`` |
| ``I_0^{*ns}`` | ``6`` | ``G_2`` | ``G_2`` | ``1`` | ``1`` | ``2`` | ``2`` | ``3`` |
| ``I_0^{*ss}`` | ``6`` | ``B_3`` | ``SO(7)`` | ``1`` | ``1`` | ``2`` | ``2`` | ``4`` |
| ``I_0^{*s}`` | ``6`` | ``D_4`` | ``SO(8)`` | ``1`` | ``1`` | ``2`` | ``2`` | ``4`` |
| ``I_1^{*ns}`` | ``7`` | ``B_4`` | ``SO(9)`` | ``1`` | ``1`` | ``2`` | ``3`` | ``4`` |
| ``I_1^{*s}`` | ``7`` | ``D_5`` | ``SO(10)`` | ``1`` | ``1`` | ``2`` | ``3`` | ``5`` |
| ``I_2^{*ns}`` | ``8`` | ``B_5`` | ``SO(11)`` | ``1`` | ``1`` | ``3`` | ``3`` | ``5`` |
| ``I_2^{*s}`` | ``8`` | ``D_6`` | ``SO(12)`` | ``1`` | ``1`` | ``3`` | ``3`` | ``5`` |
| ``I_{2k-3}^{*ns}`` | ``2k+3`` | ``B_{2k}`` | ``SO(4k+1)`` | ``1`` | ``1`` | ``k`` | ``k+1`` | ``2k`` |
| ``I_{2k-3}^{*s}`` | ``2k+3`` | ``D_{2k+1}`` | ``SO(4k+2)`` | ``1`` | ``1`` | ``k`` | ``k+1`` | ``2k+1`` |
| ``I_{2k-2}^{*ns}`` | ``2k+4`` | ``B_{2k+1}`` | ``SO(4k+3)`` | ``1`` | ``1`` | ``k+1`` | ``k+1`` | ``2k+1`` |
| ``I_{2k-2}^{*s}`` | ``2k+4`` | ``D_{2k+2}`` | ``SO(4k+4)`` | ``1`` | ``1`` | ``k+1`` | ``k+1`` | ``2k+1`` |
| ``IV^{*ns}`` | ``8`` | ``F_4`` | ``F_4`` | ``1`` | ``2`` | ``2`` | ``3`` | ``4`` |
| ``IV^{*s}`` | ``8`` | ``E_6`` | ``E_6`` | ``1`` | ``2`` | ``2`` | ``3`` | ``5`` |
| ``III^*`` | ``9`` | ``E_7`` | ``E_7`` | ``1`` | ``2`` | ``3`` | ``3`` | ``5`` |
| ``II^*`` | ``10`` | ``E_8`` | ``E_8`` | ``1`` | ``2`` | ``3`` | ``4`` | ``5`` |
| non-min. | ``12`` |  |  | ``1`` | ``2`` | ``3`` | ``4`` | ``6`` |



## Constructors

We aim to provide support for global Tate models over the following bases:
* a toric variety,
* a toric scheme,
* a (covered) scheme.
Often, one also wishes to obtain information about a global Tate model without
explicitly specifying the base space. Also for this application, we provide support.
Finally, we provide support for some standard constructions.

Before we detail these constructors, we must comment on the constructors over toric base
spaces. Namely, in order to construct a global Tate model as a hypersurface in an ambient
space, we first wish to construct the ambient space in question. For a toric base, one way to
achieve this is to first focus on the Cox ring of the toric ambient space. This ring must be
graded such that the Tate polynomial is homogeneous and cuts out a Calabi-Yau
hypersurface. Given this grading, one can perform a triangulation task. Typically, this
combinatorial task is very demanding, consumes a lot of computational power and takes a
long time to complete. Even more, it will yield a large, often huge, number of candidate
ambient spaces of which the typical user will only pick one. For instance, a common and
often appropriate choice is a toric ambient space which contains the toric base space in a
manifest way.

To circumvent this very demanding computation, our toric constructors operate in the
opposite direction. That is, they begin by extracting the rays and maximal cones of the chosen
toric base space. Subsequently, those rays and cones are extended to form one of the many
toric ambient spaces. This proves hugely superior in performance than going through the
triangulation task of enumerating all possible toric ambient spaces. One downside of this strategy
is that the so-constructed ambient space need not be smooth.

### A toric variety as base space

We require that the provided toric base space is complete. This is a technical limitation as of now.
The functionality of OSCAR only allows us to compute a section basis (or a finite subset thereof)
for complete toric varieties. In the future, this could be extended.

However, completeness is an expensive check. Therefore, we provide an optional argument which
 one can use to disable this check if desired. To this end, one passes the optional argument
 `completeness_check = false` as last argument to the constructor. The following examples
 demonstrate this:
```@docs
global_tate_model(base::AbstractNormalToricVariety; completeness_check::Bool = true)
global_tate_model(ais::Vector{T}, base::AbstractNormalToricVariety; completeness_check::Bool = true) where {T<:MPolyRingElem}
```

### A toric scheme as base space

For the same reasons as above, the toric base must be complete. Similar to toric varieties as
bases, we can use the optional argument `completeness_check = false` to switch off the
expensive completeness check. The following examples demonstrate this:
```@docs
global_tate_model(base::ToricCoveredScheme; completeness_check::Bool = true)
global_tate_model(ais::Vector{T}, base::ToricCoveredScheme; completeness_check::Bool = true) where {T<:MPolyRingElem}
```

### A (covered) scheme as base space

This functionality does not yet exist.

### Base space not specified

This method constructs a global Tate model over a base space, where
this base space is not (fully) specified. Consequently, we simply assume
that a base space exists such that the Tate sections ``a_i`` as introduced
above do exist.

For many practical applications, one wishes to assume a further factorization
of the Tate sections ``a_i``. This has the advantage that one can engineer
singularity loci or even the singularity type over a specific locus. This is
the backbone of many F-theory constructions. For example, we could consider
the factorization:
* ``a_1 = a_{10} w^0``,
* ``a_2 = a_{21} w^1``,
* ``a_3 = a_{32} w^2``,
* ``a_4 = a_{43} w^3``,
* ``a_6 = a_{65} w^5``,
In this case, it is useful to consider the polynomial ring with indeterminates
``a_{10}``, ``a_{21}``, ``a_{32}``, ``a_{43}``, ``a_{65}`` and ``w``.
In theory, one can consider these indeterminates as local coordinate of an
auxiliary base space. Indeed, for our computer implementation the polynomial
ring with these indeterminates serve as the coordinate ring of an auxiliary *toric*
base space. Despite this auxiliary base space being toric, the predictions from
such an analysis are not limited to the world of toric varieties.

For constructions along these lines, we support the following constructor:
```@docs
global_tate_model(ais::Vector{T}, auxiliary_base_ring::MPolyRing, d::Int) where {T<:MPolyRingElem}
```


### Standard constructions

We provide convenient constructions of global Tate models over
famous base spaces. Currently, we support the following:
```@docs
global_tate_model_over_projective_space(d::Int)
global_tate_model_over_hirzebruch_surface(r::Int)
global_tate_model_over_del_pezzo_surface(b::Int)
```


## Literature constructions

Certain Tate models have been studied in the physics literature over and over again.
Thereby, these constructions became famous and some were given special names. We
aim to provide support for such standard constructions. An example of such a model is
the following:
```@docs
su5_tate_model_over_arbitrary_3d_base()
```
More generally, we support literature constructions.
```@docs
literature_tate_model(; arxiv_id::String="", equ_nr::String="", description::String="")
```


## Attributes

### Basic attributes

For all global Tate models -- irrespective over whether the base is toric or not -- we support
the following attributes:
```@docs
tate_section_a1(t::GlobalTateModel)
tate_section_a2(t::GlobalTateModel)
tate_section_a3(t::GlobalTateModel)
tate_section_a4(t::GlobalTateModel)
tate_section_a6(t::GlobalTateModel)
tate_polynomial(t::GlobalTateModel)
```
In case the global Tate model is constructed over a not fully specified base,
recall that we construct an auxiliary (toric) base space as well as an
auxiliary (toric) ambient space. The (auxiliary) base and ambient space can
be accessed with the following functions:
```@docs
base_space(t::GlobalTateModel)
ambient_space(t::GlobalTateModel)
fiber_ambient_space(t::GlobalTateModel)
```
The following method allows to tell if the base/ambient space is auxiliary or not:
```@docs
base_fully_specified(t::GlobalTateModel)
```
The user can decide to get an information whenever an auxiliary base space,
auxiliary ambient space or auxiliary hypersurface have been computed.
To this end, one invokes `set_verbosity_level(:GlobalTateModel, 1)`.
More background information is available
[here](http://www.thofma.com/Hecke.jl/dev/features/macros/).

### Advanced attributes

The following attributes are currently only supported in a toric setting:
```@docs
calabi_yau_hypersurface(t::GlobalTateModel)
weierstrass_model(t::GlobalTateModel)
```
Note that for applications in F-theory, *singular* elliptic fibrations are key
(cf. [Wei18](@cite) and references therein). Consequently the discriminant
locus as well as the singular loci of the fibration in question are of ample
importance:
```@docs
discriminant(t::GlobalTateModel)
singular_loci(t::GlobalTateModel)
```

### Attributes for literature models

For literature models, we provide the following attributes:
```@docs
doi(t::GlobalTateModel)
arxiv_id(t::GlobalTateModel)
version(t::GlobalTateModel)
equation_number(t::GlobalTateModel)
description(t::GlobalTateModel)
link(t::GlobalTateModel)
```
Certainly, one can add this information for a model that does
not have it:
```@docs
set_description(t::GlobalTateModel, description::String)
```
Note however, that these changes will (currently) not be stored
in our data base. One can also check if a model has a particular
set of information. This is achieved with the following methods:
* `has_doi(t::GlobalTateModel)`,
* `has_arxiv_id(t::GlobalTateModel)`,
* `has_version(t::GlobalTateModel)`,
* `has_equation_number(t::GlobalTateModel)`,
* `has_description(t::GlobalTateModel)`,
* `has_link(t::GlobalTateModel)`.


## Methods

### Fiber study

In F-theory, it is standard to not work with the singular space directly.
Rather, one resolves its singularities in order to obtain a smooth space
instead. Subsequently, one performs computations on this smooth space.

In order to perform such a resolution, one wishes to analyze the fibration
in detail. The following method aims at giving a first window into this analysis
by working out the fiber components and their intersection pattern over a
particular locus of the base.
```@docs
analyze_fibers(model::GlobalTateModel, centers::Vector{<:Vector{<:Integer}})
```

### Resolution(s) of a singular model

A central task in F-theory is to resolve a singular global Tate model.
For literature models, we have stored resolutions in our data base.
Upon construction of a literature model, we load these known resolutions.
The user can list all the resolutions in our database as follows:
```@docs
resolutions(t::GlobalTateModel)
```
In addition, some user might want to add a resolution. This can be achieved with
the following method:
```@docs
add_resolution(t::GlobalTateModel, centers::Vector{Vector{String}}, exceptionals::Vector{String})
```
Provided that a resolution for a model is known, we can (attempt to) resolve the model.
```@docs
resolve(t::GlobalTateModel, index::Int)
```
