```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["space_germs.md"]
```

# Space Germs

## [Generalities on Space germs](@id space_germ_generalities)
 
The geometric notion of a space germ is a local concept. A space germ $(X,x)$ at a point $x$ is an equivalence class of ringed spaces, each of which contains $x$ in its underlying topological space, and the equivalence relation is precisely the existence of an open neighbourhood of $x$ on which the spaces coincide.

Depending on the kind of ringed space in question, space germs arise in
different forms:

  * a space germ in the context of affine schemes is the geometric object arising from a given scheme by localization at a point, leading to the stalk of the structure sheaf at the respective prime ideal.

  * in the context of singularity theory, the (anti-)equivalence of categories between complex space germs and analytic ``{\mathbb CC}``-algebras allows the direct definition of a space germ from algebraic data

Note that analytic algebras as mentioned above, have two computational problems.
On one hand, exact computations can only be performed over fields in OSCAR
permitting exact computations, in particular not over ``{\mathbb R}`` or 
``{\mathbb C}``. This usually does not pose a problem, if the input data is 
in an exact smaller field. But unfortunately, also analytic algebras themselves 
do not allow exact computations so that applications have to be considered 
in a localization of an affine algebra. With due care, output may again be 
interpreted in terms of a multivariate formal power series using the following 
inclusions:
```math
{\mathbb Q}[\underline{x}]_{\langle x \rangle} \hookrightarrow
  {\mathbb Q}\{\underline{x}\}\hookrightarrow
  {\mathbb Q}[[\underline{x}]]
```
At particular points, where this difficulty of interpretation 
manifests itself prominently, a suitable warning, note or example has been 
placed (but certainly not everywhere). 

Textbooks covering space germs in the sense of analytic algebras and
singularity theory are:
- [GLS07](@cite)
- [JP00](@cite)

For the point of view of schemes, we refer to the page on schemes and the references given there; an example of a standard textbook is
- [Har77](@cite).

## Creating Space Germs in OSCAR

In general, space germs in OSCAR are created in the following ways:

 * localization of an affine scheme at a point

   ```julia
   SpaceGerm(X::AbsSpec, I::Ideal)
   ```
   where I is a (maximal) ideal describing the chosen ${\mathbb k}$-point on the affine ${\mathbb k}$-scheme X.  
   **Provides:** SpaceGerm.
   ```julia
   germ_at_point(X::AbsSpec, I::Ideal)
   ```
   where I is a (maximal) ideal describing the chosen ${\mathbb k}$-point on the affine ${\mathbb k}$-scheme.  

   **Provides:** SpaceGerm, restriction map.


 * localization of a polynomial ring at a point

   ```julia
   germ_at_point(R::MPolyRing, I::MPolyIdeal)
   germ_at_point(R::MPolyQuo, I::MPolyIdeal)
   ```
   where I is a maximal ideal describing the chosen base_ring(R)-point.

   **Provides:** space germ and restriction map.

 * localized ring with respect to complement of prime ideal or complement of maximal ideal)

   ```julia
   SpaceGerm(R::MPolyLocalizedRing)
   SpaceGerm(R::MPolyQuoLocalizedRing)
   ```
   **Provides:** Space germ
   ```julia
   germ_at_point(R::MPolyLocalizedRing)
   germ_at_point(R::MPolyQuoLocalizedRing)
   ```

   **Provides:** Space germ, restriction map  

   (point inherited from underlying local ring)

As for all affine schemes, morphisms of space germs are in $1:1$ correspondence with morphisms of the underlying local rings. Hence
also handled in this way in OSCAR.

## Basic functionality for space germs

Most of the basic functionality immediately falls back to the underlying 
affine algebra or affine scheme and is just provided on the geometric side 
for convenience and consistence of functionality:

### internal data of a space germ  

 * Pass from the germ $(X,x)$ back to some affine scheme $X=Spec R$, where $R$ is a quotient of a multivariate polynomial ring and localizes to ${\mathcal O}_{X,x}$.
   ```julia
   representative(X::SpaceGerm)
   ```
   **Provides:** AffineScheme   

!!! note
    At first need or at the latest upon the first explicit call of representative, the respective affine scheme is cached and subsequently used for all further purposes requiring a representative.

 * Given the space germ $(X,x)$, the point $x$ is returned by:

   ```julia
   point(X:SpaceGerm)
   ```
   **Provides:** Vector describing of point coordinates
 
 * Given a space germ $(X,x)$, the corresponding local ring ${\mathcal O}_{(X,x)}$ is returned by:
   ```julia
   ring(X::SpaceGerm)
   ```
   **Provides:** `MPolyQuoLocalizedRing` or `MPolyLocalizedRing`

 * Analogously, the modulus of the ring of a given germ $(X,x)$ can be obtained by:
   ```julia
   ideal(X::SpaceGerm)
   ```
   **Provides:** Ideal

 * For technical reasons all germs $(X,x)$ in OSCAR are embedded and the smooth ambient germ can be accessed by:
   ```julia
   ambient_germ(X::SpaceGerm)
   ```
   **Provides:** SpaceGerm

### containment/equality of space germs
    
 * containment  

   ```julia
   issubset(X::SpaceGerm, Y::SpaceGerm)
   ```

   Test whether X is a subgerm of Y.  

   **Provides:** Boolean Value
 
!!! note
    The name 'issubset' has been chosen for consistency with other types, but is slightly misleading, as the function does not only test containment of the underlying sets, but in the scheme-theoretic sense.

 * equality

   ```julia
   X ==Y
   ```

   **Provides:** Boolean Value

!!! note
    Equality of the germs X and Y is tested in the sense of equality of subgerms of the same ambient space germ.  


 * emptyness

   For the test of equality to the empty germ, a separate function is available:
   ```julia
   is_empty(X::SpaceGerm)
   ```

   **Provides:** Boolean Value

 * intersection

   Compute the intersection of two subgerms of a common ambient germ
   ```julia
   intersect(X::SpaceGerm,Y::SpaceGerm)
   ```
     **Provides**: SpaceGerm   

### singular locus

```julia
singular_locus(X::SpaceGerm)
```

   Computes the germ of the singular locus of the given germ.

   **Provides:** SpaceGerm, inclusion morphism  
  
If a test for smoothness is the goal, the following functions can be used

```julia
is_smooth(X::SpaceGerm)
```

   **Provides:** Boolean Value

<!-- ```julia is_regular(X::SpaceGerm)```  still needs to be implemented -->

## Identification of important classes of space germs

The following functions set/read properties which identify particular cases. 
Certain properties and invariants are only defined for one or several special
classes, whereas some other properties and invariants make sense for arbitrary
space germs, but may only be accesible to (efficient) computations in 
particular classes. The information on the currently supported particular cases
* [Hypersurface Singularities](@ref hypersurface_sing)
* [Complete Intersections Singularities](@ref ICIS_sing)
* Cohen-Macaulay Codimension 2 Singularities 
can be found on the respective pages.

```julia
is_isolated(X::SpaceGerm)
is_hypersurface(X::SpaceGerm)
is_complete-intersection(X::SpaceGerm)
is_CMC2(X::SpaceGerm)
```

Additionally, rigid singularities are detectable:
``` julia
is_rigid(X::SpaceGerm)

```

**Provides:** Boolean Value

!!! note  
    The properties isolated hypersurface singularity (IHS), isolated complete intersection singularity (ICIS) and isolated Cohen-Macaulay codimension 2 singularity (ICMC2) may each be deduced from the appropriate combinations of the aforelisted functions.

!!! note "on planned implementation"
    Use minimal set of generators for hypersurface and complete intersection, use Hilbert-Burch structure for CMC2 -- store minimal set of generators or presentation matrix. For rigidity test use T1_module below

## Further basic functionality for arbitrary space germs

### Hilbert-Samuel polynomial

The Hilbert-Samuel polynomial w.r.t. the maximal ideal
```julia
Hilbert-Samuel-polynomial(X::SpaceGerm)
```

**Provides:** MPolyElem

### Dimension

```julia
dimension(X::SpaceGerm)
```
**Provides:** integer value

### multiplicity

```julia
multiplicity(X::SpaceGerm)
```
**Provides**: integer value

### associated graded ring

```julia
associated_graded_ring(X::SpaceGerm)
```
Determines the associated graded ring to the local ring of the space germ w.r.t. its maximal ideal.

**Provides:** graded ring

!!! note
    'dimension' and 'multiplicity' both rely on data from the Hilbert-Samuel polynomial.

### Decomposition into irreducible components

Based on primary [decomposition of ideals](@ref ideal_decomp) in multivariate
polynomial rings, a decomposition of a space germ into irreducible components 
may be computed (in favorable settings). We only mention this topic here to 
point out the pitfalls in interpreting the result in terms of complex space germs. These problems are only due to the fact that we work over computable OSCAR fields and in localizations of polynomial rings.  

**Pitfall** **1:**  
$x^2+y^2 = (x+iy)(x-iy)$ is obviously reducible over the complex numbers ir even ${\mathbb Q}[i]$, but irreducible over the rationals.
As outlined in [decomposition of ideals](@ref ideal_decomp) and explained in the literature cited there, this can be detected/avoided by using absolute primary decomposition.

**Pitfall** **2:**  
Even in cases, where pitfall 1 does not cause any problems, we might still not be able to separate all components of a space germ computationally. This is due to the fact that the inclusion of rings 
${\mathbb Q}[\underline{x}]_{\langle x \rangle} \subsetneq {\mathbb Q}\{\underline{x}\}$ is strict and computations are only possible in the smaller one, 
whereas results often need to be interpreted in the larger one. For instance,
the polynomial $y^2-x^2-x^3$ is irreducible in the smaller ring, but factors as
$(y-x(\sqrt{1+x})\cdot(y+x(\sqrt{1+x})$ in the larger one.

## [Deformations of space germs](@id spacegerms_deform)

### Preliminary remarks on deformations

To provide the context and, in particular, fix notation for the deformation theoretic functionality which follows, we recall some notions and notation from [GLS07](@cite), Chapter II. (Note that the use of certain terms differs between textbooks in this field.)  

Let $(X,x)$ and $(S,s)$ be space germs. A deformation of $(X,x)$ over $(S,s)$ is a flat morphism $\phi$ of germs together with an isomorphism from $(X,x)$ to the fibre $({\mathcal X}_s,{\mathcal x})$ of $\phi$, best visualized by a cartesian diagram
```math
(X,x) \stackrel{\iota}{\hookrightarrow} ({\mathcal X},{\mathcal x})\\
\phantom{(X,)}\downarrow \phantom{\;\;\;\longrightarrow\;\;\;} \;\; \downarrow \phi \quad \textrm{flat}\\
\{pt\} \quad \hookrightarrow (S,s)
```
For convenience of notation, we will refer to this deformation as the deformation $(\iota,\phi)$, as all information is already encoded in these two morphisms.

Given two deformations $(\iota_1, \phi_1)$ and $(\iota_2,\phi_2)$ of $(X,x)$ over $(S,s)$, a morphism of deformations from $(\iota_1,\phi_1)$ to $(\iota_2,\phi_2)$  consists of a pair of morphisms $(\psi,\varphi)$, where $\psi: ({\mathcal X}_1,{\mathcal x}_1) \longrightarrow ({\mathcal X}_2,{\mathcal x}_2)$ is a morphism of the total spaces and $\varphi: (S_1,s_1) \longrightarrow (S_2,s_2)$ a morphism of the base spaces, such that the diagram below commutes:

BILD EINFUEGEN -- DIAGRAMM HAT SCHRAEGE PFEILE...  

Given a deformation $(\iota,\phi)$ of $(X,x)$ over $(S,s)$ and a morphism of
germs $\varphi: (T,t) \longrightarrow (S,s)$, then the base change map 
$\varphi$ induces a deformation $\varphi^*(\iota,\phi)=(\varphi^*\iota,\varphi^*\phi)$ of $(X,x)$ over $(T,t)$ in the following way:

BILD EINFUEGEN -- DIAGRAMM HAT SCHRAEGE PFEILE...

### Versality

A deformation $(\iota,\phi)$ of $(X,x)$ over $(S,s)$ is called complete, if 
any other deformation $(\iota_2,\phi_2)$ of $(X,x)$ over some $(T,t)$ arises 
from it by suitable base change. A complete deformation $(\iota,\phi)$ of 
$(X,x)$ is called (formally) versal, if for any deformation $(\iota_2,\phi_2)$ 
of $(X,x)$ over some base space $(T_2,t_2)$, any closed embedding of 
complex (artinian) germs $\kappa: (T_1,t_1) \hookrightarrow (T_2,t_2)$ and 
any morphism $\varphi_1:(T_1,t_1) \longrightarrow (S,s)$ satisfying 
$(\varphi_1^*\iota,\varphi_1^*\phi) \cong (\kappa^*\iota_2,\kappa^*\phi_2)$ 
there exists a morphism $\varphi: (T_2,t_2) \longrightarrow (S,s)$ such that 
the following following diagram of Cartesian squares holds:


BILD EINFUEGEN -- DIAGRAMM HAT SCHRAEGE PFEILE...  

A (formally) versal deformation is called semiuniversal or miniversal, if the
Zariski tangent map $T(\varphi): T_{(T_2,t_2)} \longrightarrow T_{(S,s)}$ is
uniquely determined by $(\iota,\phi)$ and $(\iota_2,\phi_2)$. By a famous 
result of Grauert, any complex space germ with an isolated singularity 
possesses a semiuniversal deformation.

 Moreover, the main intermediate results, i.e. $T^1_{(X,x)}$ and $T^2_{(X,x)}$, are of interest of their own and will be discuss below.

``` julia
versal_deformation(X::SpaceGerm)
versal_deformation(X::SpaceGerm,n::Integer)
```

**Provides:** pair of morphims

!!! note
    For isolated [hypersurface singularities](@hypersurface_deform) and isolated complete intersection singularities (HIER LINK SETZEN SOBALD SEITE EXISTIERT) the computation of versal deformations is straight forward and described in the respective sections.  In the general case, the computation is based on an iteration, whose termination is not guaranteed. The standard precision of the computed approximate result (set to ???) may be overwritten by the additional argument. 

To help understand the meaning of the above mentioned precision $n$, we very briefly sketch the key idea of the iteration: As a starting point, deformations over the smallest fat point, corresponding to ${\mathbb C}[\varepsilon]/\langle \varepsilon^2\rangle$, are considered (cf. also $T^1_{(X,0)}$ below) and then iteratively lifted to larger and larger fat points as base space. In good situations (such as hypersurfaces and ICIS) this iteration stops after finitely many steps; in unlucky situations it is indeed an infinite process. Hence this provides a formally versal deformation over a base space with local ring ${\mathbb C}[[\underline{t}]]/J$ for some ideal $J$, but in practice computations have to stop at some finite order to finish in finite time. 

### First order deformations

The module of first order deformations can be understood as the Zariski tangent space to the base of the versal deformation.

Denote by $(T_{\varepsilon},0)$ the complex space germ consisting of a point and a local ring ${\mathbb C}[\varepsilon]/\langle \varepsilon^2\rangle$. For a space germ $(X,0) \subset ({\mathbb C}^n,0)$ the ${\mathbb C}$-vector space of isomorphism classes of deformations of $(X,0)$ over $T_{\varepsilon}$ is denoted by $T^1_{(X,0)}$ and its elements are usually referred to as first order infinitessimal deformations of $(X,0)$. Using the notation ${\mathcal O}_{X,0} = {\mathcal O}_{{\mathbb C}^n,0}/I$, it can be computed from the normal module ${\mathcal N}_{X/{\mathbb C}^n,0} \cong Hom_{{\mathcal O}_{{\mathbb C}^n,0}}(I,{\mathcal O}_{{\mathbb C}^n,0})$ via the following exact sequence:
```math
0 \longrightarrow \theta_{X,0}
\longrightarrow \theta_{{\mathbb C}^n,0} \otimes_{{\mathcal O}_{{\mathbb C}^n,0}} {\mathcal O}_{X,0}
\longrightarrow {\mathcal N}_{X/{\mathbb C}^n,0}
\longrightarrow T^1_{(X,0)} \longrightarrow 0
```
where $\theta_{X,0}$ denotes the sheaf of ${\mathbb C}$-derivations of ${\mathcal O}_X$ with values in  ${\mathcal O}_X$ and $\theta_{{\mathbb C}^n,0}$ the one of ${\mathcal O}_{{\mathbb C}^n}$ with values in ${\mathcal O}_{{\mathbb C}^n}$, each locally at $0$.

``` julia
T1_module(X::SpaceGerm)
```

**Provides:** module

``` julia
T1_basis(X::SpaceGerm)
```

**Provides:** a pair consisting of a vector and a list of vectors

!!! note
    The desired object $T^1_{(X,0)}$ carries the structure of an an ${\mathcal O}_{(X,0)}$-module and of an $O_{({\mathbb C}^n,0)}$-module. It is returned in the latter way by T1-module. It also carries the structure of a ${\mathbb C}$-vector space, which may or may not be finite-dimensional for a given singularity. If its dimension is not finite, T1-basis returns an error. 

If the vector space dimension of $T^1_{(X,0)}$ is finite, T1-basis returns the following data:
* entry 1: $\underline{f} = (f_1,\dots,f_k) \in {\mathcal O}_{({\mathbb C}^n,0)}^k$ such that ${\mathcal O}_{(X,0)} =  {\mathcal O}_{({\mathbb C}^n,0)}/ \langle f_1,\dots,f_k \rangle$ 
* entry 2: vectors $g_1,\dots,g_{\tau}$ in ${\mathcal O}_{({\mathbb C}^n,0)}^k$, each with a single monomial entry such that 
``` math
T^1_{(X,0)} \cong \bigoplus_{i=1}^{\tau} {\mathbb C} g_i 
```
via the identification of an element $\sum_{i=1}^{\tau} a_i g_i$ of the right hand side with the isomorphism class of the deformation $(X,0) \hookrightarrow ({\mathcal X},0) \longrightarrow (T_{\varepsilon},0)$ defined by $\langle F_1,\dots,F_k \rangle$ where


``` math
\underline{F} = (F_1,\dots,F_k) = \underline{f} + \varepsilon \sum_{i=1}^{\tau} a_ig_i
```
(see [JP00](@cite) for a detailed discussion including an example).

!!! note "on planned implementation"
    T1_module mimics the steps from sing.lib, T1-basis uses T1_module

### Obstructions

The flatness condition on the map $\phi$ in a deformation $(\iota, \phi)$ implies that lifting a deformation over $(T
_{\varepsilon})$ to a deformation over a higher order fat point is not always feasable. The module $T^2_{(X,0)}$ desc
ribes the obstructions for lifting a deformation over $(T_{\varepsilon},0)$ to one over an infinitessimally bigger on
e $(T',0)$.

More precisely, $T^2_{(X,0)}$ is defined in the following way:
Let
``` math
{\mathcal O}_{({\mathbb C}^n,0)}^\ell \longrightarrow 
{\mathcal O}_{({\mathbb C}^n,0)}^k \stackrel{\varphi_f}{\longrightarrow} 
\langle I_{(X,0)} \rangle
```
be a free presentation of $I_{(X,0)}=\langle f_1,\dots,f_k \rangle$. Denote by $M_R \subset {\mathcal O}_{({\mathbb C}^n,0)}^k$ the kernel of $\varphi_f$ and by $M_K \subset {\mathcal O}_{({\mathbb C}^n,0)}^k$ the module of Koszul-relations among the $f_1,\dots,f_k$. Then the inclusion $M_R \subset {\mathcal O}_{({\mathbb C}^n,0)}^k$ and the fact that all Koszul-relations are zero modulo $I_{(X,0)}$ provide a ${\mathcal O}_{(X,0)}$-linear map $M_R/M_K \longrightarrow {\mathcal O}_{(X,0)}^k$. Taking the dual we obtain the exact sequence which defines $T^2_{(X,0)}$:
``` math
Hom_{{\mathcal O}_{(X,0)}}({\mathcal O}_{(X,0)}^k, {\mathcal O}_{(X,0)})
\longrightarrow
Hom_{{\mathcal O}_{(X,0)}}(M_R/M_K,  {\mathcal O}_{(X,0)})
\longrightarrow T^2_{(X,0)} \longrightarrow 0
```

It can be computed using

``` julia
T2_module(X::SpaceGerm)
```

**Provides**: module

``` julia
T1T2_modules(X::SpaceGerm)
```

**Provides**: pair of modules


!!! note
    If $T^2_{(X,0)} = 0$, like for hypersurfaces and complete intersections, any first order deformation can be lifted and the base of the semiuniversal deforation is smooth. An explicitly computed example with non-smooth base of the versal deformation can be found in [JP00](@cite), Chapter 10.2.

!!! note
    If $T^1$ and $T^2$ are both needed, it is more efficient to compute them simultaneously via the function T1T2_modules. 

### Families of space germs

Given a pair of maps $(\iota,\phi)$ satisfying a commutative diagram as above, flatness of $\phi$ can be checked explicitly via the relation-lifting-property. This is available as

``` julia
is_flat_deformation(iota::MapGerm,phi::MapGerm)
```

**Provides:** Boolean value

For a given deformation $(\iota,\phi)$ versality and minimality can be checked via the Kodaira-Spencer map. This is available as

``` julia
is_versal(iota::MapGerm,phi::MapGerm)
is_semiuniversal(iota::MapGerm,phi::MapGerm)
```

**Provides:** Boolean value

!!! note
    Recognition of trivial subfamilies via the Kodaira-Spencer map is available for hypersurfaces, ICIS and ICMC2 singularities and described in the respective sections.

