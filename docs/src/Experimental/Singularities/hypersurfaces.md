```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["hypersurfaces.md"]
```

# Hypersurface Singularities

Hypersurface singularities are space germs of the form $(V(f),p)$ for some power series f at a point $p$. As explaind in section [Generalities on Space Germs](@ref space_germ_generalities), the need for exact computations forces OSCAR to restrict all considerations to the subring of those power series representable by means of an element of the localization of an affine scheme at a point. 

Where no other reference is specified, we refer the reader to the textbook [GLS07](@cite).

For simplicity of notation, ${\mathcal O}_{n,p}$ denotes the localization of the polynomial ring in $n$ variables at the point $p$ here and ${\mathfrak m}_p$ denotes its maximal ideal. If $p$ is origin, we often omit $p$ in the index.

**Right** **equivalence** of hypersurface singularities takes the point of 
view of the a map $f: ({\mathbb C}^n,p) \mapsto ({\mathbb C},0)$ defined by the power series $f(x)$ 
and allows coordinate changes in the source,i.e.
```math
f \sim_R g :\Longleftrightarrow \;\exists \varphi \in Aut({\mathcal O}_{n,p}): f = g \circ \varphi.
```
**Contact** **equivalence** takes the point of view of the vanishing locus of f and
allows isomorphism of space germs, i.e. 
```math
f \sim_C g :\Longleftrightarrow (V(f),p) \cong (V(g),p).
```
 
Due to these two perspectives, the majority of functionality for hypersurfaces is available both for elements of ${\mathcal O}_{n,0}$ and for space germs at an arbitrary point $p \in {\mathbb C}^n$ with coordinates in a computable subfield $k$ of ${\mathbb C}$. To keep the notation simple and cover both cases simultaneously, the following discussion will be focusing on germs at the origin. However, the functions accepting a space germ as input may be called with a hypersurface germ at any (k-)point $p$.

!!! note
    All functionality for space germs is available for hypersurface germs, but several algorithms possess significantly more efficient variants exploiting the knowledge of handling a hypersurface germ.

## Basic properties of Hypersurface Singularities

### isolatedness and singular locus

A hypersurface singularity (V(f),0) is called isolated, if it posseses a 
representative on a sufficiently small neighbourhood of 0, whose singular
locus is the origin. By abuse of notation, we also refer to the germ of
this singular locus as the singular locus of the original germ.

```julia
   is_isolated(f::MPolyElemLoc)
   is_isolated(X::HypersurfaceGerm)
```

**Provides:** Boolean value

```julia
   singular_locus(f::MPolyElemLoc)
   singular_locus(X::HypersurfaceGerm)
```

**Provides:** SpaceGerm

!!! note "on planned implementation"
    falls back to space germ, using its singular_locus and is_isolated

### Finite determinacy

A hypersurface germ $(V(f),0)$ is called finitely R-determined (or finitely
C-determined), if 
```math
f \equiv g \;\; mod \;\; m^{k+1} \Longrightarrow f \sim_R g \;\;(or \;\; f \sim_C g \;\;respectively).
```

```julia
is_finitelydetermined(f::MPolyElemLoc)
is_finitelydetermined(X::HypersurfaceGerm)
```

**Provides:** boolean value


```julia
determinacy_bound(f::MPolyElemLoc)
determinacy_bound(X::HypersurfaceGerm)
```

**Provides:** integer

Return some (in general non-optimal) determinacy bound, if the input is finitely determined; otherwise, it returns the negative integer -1.

**Provides:** integer

!!! note "on planned implementation"
    use highest corner of Milnor algebra. 

### Quasihomogeneity

A hypersurface germ $(X,0)$ is called quasihomogeneous, if it possesses a quasihomogeneous representative $V(f)$. More precisely, there are an integer $d$ and a tuple of positive rational numbers $w = (w_1,\dots,w_n)$ such that all monomials in $f$ are of $w$-weighted degree $d$, i.e. $\sum_{i=1}^n w_ia_i =d$  for any monomial $x^a=x_1^{a_1}\cdots x_n^{a^n}$ in the support of $f$. We will refer to these data as a weight system $(d,w)$.

```julia
is_quasihomogeneous(f::MPolyElemLoc)
is_quasihomogeneous(X::HypersurfaceGerm)
```

**Provides:** boolean value

quasihomogeneous_weights(f::MPolyElemLoc)
quasihomogeneous_weights(X::HypersurfaceGerm)

**Provides:** pair of integer and vector of rationals

If the germ is quasihomogeneous w.r.t. a weight system $(d,w)$, this pair is
returned. Otherwise, the pair $(0,[0,\dots,0])$ is returned.

### Morse critical point

A germ $(V(f),0)$ is called a Morse critical point, if it is singular and
the Hessian Matrix $H(f)$ of $f$ at $0$ has full rank.

```julia
   is_morse_critical_point(f::MPolyElemLoc)
   is_Morse_critical_point(X::HypersurfaceGerm)
```

**Provides:** Boolean Value

## Newton polytope and related data

!!! note  
    The functionality for Newton polytopes is only provided for f and not for hypersurface germs to avoid confusion, because the construction depends on the choice of $f$ or of the embedding of the germ into the ambient germ.

!!! note "on planned implementation"
Implementation is already underway as sideeffect of some thesis. Do not touch at the moment.

### Newton polytope

For a multivariate power series $f = \sum_{a \in NN^n}$ c_a x^a$, the Newton
polytope is defined as

```math
\Delta (f) = Conv \{a \in NN^n \mid c_a \neq 0\}
```

```julia
   Newton_polytope(f::MPolyElemLoc)
```

**Provides:**              

### Newton diagram

Given the Newton polytope $\Delta(f)$ of a power series $f$, the Newton diagram
$\Gamma(f)$ arises from the following construction:
```math
K(f)=Conv( \{ 0 \} \cup \Delta(f))\\
K_0(f)=\overline{K(f) \setminus \Delta(f)}\\
\Gamma(f)= K_0(f) \cap \Delta(f)

```
The Newton diagram can also be understood as the union of the bounded faces of 
```math
Conv(\bigcup_{a \in supp(f)} a+{\mathbb R}_{\geq 0}^n)
```

!!! note
    This command is only provided for $f$ and not for space germs as the 
    construction depends on the chosen system of parameters.


```julia
   Newton_diagram(f::MPolyElemLoc)
```

**Provides:** 

### is_convenient

A power series f is called convenient, if its Newton Diagram meets all the
coordinate axes. The implemented algorithm does not explicitly compute the 
Newton diagram!

```julia
   is_convenient(f::MPolyElemLoc)
```

**Provides:** Boolean value

!!! note "on planned implementation "
    This does no require computing the Newton-diagram. It is a matter of evaluating f on the axes.

### is_Newton-non-degenerate

Let $f = \sum_{a \in NN^n}$ c_a x^a$ be a multivariate, convenient power series and let $\Gamma(f)$ be its Newton-diagram. For each face $\sigma \in \Gamma(f)$ we define a quasihomogeneous polynomial 
```math
   f_{\sigma}=\sum_{a \in \sigma} c_a x^a
```
f is called Newton-non-degenerate, if for all faces $f_{\sigma}$, $\sigma \in \Gamma(f)$, the germ $(V(f_\sigma),0)$ is non-singular outside the coordinate hyperplanes $(V(\prod_{i=1}^n x_i),0)$.

```julia
   is_Newton_non_degenerate(f::MpolyElemLoc)
```

**Provides:** Boolean Value

## Invariants of Hypersurface Singularities

Several geometrical, topological and deformation theoretical properties of hypersurface singularities can be described by numerical invariants. Those which are accessible to computation within OSCAR are listed here

### corank

For a singular hypersurface germ $(V(f),0) \subset ({\mathbb C}^n,0)$, the corank is the integer $n-rank(H(f))$ where $H(f)$ denotes the Hessian matrix of $f$, i.e. the $n \times n$ matrix of second order partial derivatives.  

By the Generalized Morse Lemma the following holds for a space
germ $(V(f),0) \subset ({\mathbb C}^n,0)$ of corank $c$:
```math
   \exists g \in \langle x_1,\dots,x_c \rangle^3 \subset {\mathbb C}\{x_1,\dots,x_c\} : f \sim_R g(x_1,\dots,x_c) + \sum_{i=c+1}^n x_n^2
```

```julia
   corank(f::MPolyElemLoc)
   corank(X::HypersurfaceGerm)
```

**Provides:** integer

!!! note "on planned implementation"
    this is 3 lines, implement in OSCAR directly

### multiplicity

For a hypersurface germ $(V(f),0)$ the multiplicity at the origin is the
order of the power series expansion of $f$, i.e.
```math
   mult(V(f),0) = ord_0(f) = min\{k\in {\mathbb N} \mid f \not\in m^{k+1}\}
```

```julia
   multiplicity(f::MPolyElemLoc)
   multiplicity(X::HypersurfaceGerm)
``` 

**Provides:** integer

### Milnor number and Milnor Algebra

For an isolated singularity $(V(f),0) \subset ({\mathbb C},0)$, its Milnor fibre has the topological type of a bouquet of $(n-1)$-spheres and there are precisely $\mu$ such spheres. This number $\mu$ is called the Milnor number of $(V(f),0)$ and can be computed as the vector space dimension of the Milnor algebra:  

```math
\mu(V(f),0) = dim_{{\mathbb C}} ({\mathcal O}_n / \langle \frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

**Disambiguation:**
For an affine hypersurface $V(f)$ possessing at most isolated singularities, 
the sum over the Milnor numbers at all singular points is called its global Milnor number and can be computed as  

```math
\mu(V(f)) = dim_{{\mathbb C}} ({\mathbb C}[x_1,\dots,x_n]  / 
             \langle \frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

```julia
   milnor_number(f::MPolyElemLoc)
   milnor_number(X::HypersurfaceGerm)
   milnor_number(f::MPolyElem)
   milnor_number(X::AffineScheme)
```

**Provides:** integer


```julia
   milnor_algebra(f::MPolyElemLoc)
   milnor_algebra(X::HypersurfaceGerm)
   milnor_algebra(f::MPolyElem)
   milnor_algebra(X::AffineScheme)
``` 

**Provides:** MPolyQuo or MPolyQuoLocalizedRing (depending on input)

!!! warning  
    In contrast to other functions for space germs, milnor_number and 
    milnor_algebra will not throw an error when called with an affine 
    scheme or a polynomal ring element, but compute a different value!

!!! note "on planned implementation"
    implement in OSCAR, good simple use case for interplay between 0-dim rings and fin-dim CC-vector spaces.

### Tjurina number

For an isolated singularity $(V(f),0) \subset ({{\mathbb C}}^n,0)$, the 
base of a miniversal family has dimension $\tau$. This number $\tau$ is 
called the Tjurina number of $(V(f),0)$ and can be computed
as the vector space dimension of the Tjurina algebra:

```math
   \tau(V(f),0) = dim_{{\mathbb C}} ({\mathcal O}_n / \langle f,\frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

For an affine hypersurface $V(f)$ possessing at most isolated singularities, 
the sum over the Tjurina numbers at all singular points is called the global Tjurina number and can be computed as

```math
   \tau(V(f)) = dim_{{\mathbb C}} ({{\mathbb C}}[x_1,\dots,x_n]  / 
             \langle f,\frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

```julia
   tjurina_number(f::MPolyElemLoc)
   tjurina_number(X::HypersurfaceGerm)
   tjurina_number(f::MPolyElem)
   tjurina_number(X::AffineScheme)
``` 

**Provides:** integer

!!! warning  
    In contrast to other functions for space germs, tjurina_number will not throw an error when called with an affine scheme or a polynomal ring element, but compute a different value.

!!! note
    Please also refer to the section on $T^1$ in the section space germs for a different perspective on the Tjurina algebra.

!!! note "on planned implementation"
    implement in Oscar , same comment as for Milnor number -- should be implemented in analogous way. 
 
### monodromy

For a detailed introduction to the monodromy of a hypersurface singularity we
refer the reader to [C-MTS20](@cite), chapter 8. Here we assume familiarity 
with the notion and explain what data is provided by the respective functions
in OSCAR.  

Denote by $X_t$ a Milnor fibre of a $n$-dimensional isolated hypersurface 
singularity $(X,0) \subset ({\mathbb C}^{n+1},0)$. Then the Monodromy Theorem states
that the eigenvalues of the monodromy $M \in Aut(H^n(X_t),ZZ)$ are roots of
unity and its Jordan blocks are of size at most $(n+1) \times (x+1)$. The
eigenvalues and the sizes and multiplicities of the corresponding Jordan 
blocks can be determined using OSCAR:

```julia
   monodromy(f::MPolyElemLoc)
   monodromy(X::HypersurfaceGerm)
```

**Provides:** list of tuples, each of the form (eigenvalue, size of block, multiplicity of block)

!!! note  
    The eigenvalues are not returned as roots of unity, but as a rational
    number $q \in {\mathbb Q}$ such that $e^{2 \pi i q}$ is the desired eigenvalue. 

!!! note "on planned implementation"
    wrap gmssing.lib. mondromy.lib is older, less general and in some cases less stable (communicate with Matthias Schulze before beginning)


### spectrum and spectral pairs

For a detailed treatment of the spectrum of a singularity and of the spectral pairs, we refer the reader to the textbook [AGV85](@cite), volume II. It is a powerful invariant, which permits to exclude equivalence or adjacency of given hypersurface germs in many practical cases.  

In OSCAR the singularity spectrum can be computed using two different methods:

```julia
   spectrum(f::MPolyElemLoc)
   spectrum(X::HypersurfaceGerm)
   spectrumNND(f:MPolyElemLoc)
```

 **Provides:** list of pairs (rational, integer) encoding the spectral numbers and their multiplicities  

!!! note  
    spectrumNND is only applicable, if the input $f$ is Newton non-degenerate. This is checked at the start of the procedure and an error is returned, if $f$ does not have this property. If spectrumNND is applicable, it uses a combinatorial algorithm and is often faster than the function spectrum which relies on the Brieskorn lattice/the Gauss-Manin connection.

```julia
   spectral-pairs(f::MpolyElemLoc)
   spectral-pairs(X::HypersurfaceGerm)
```

**Provides:** list of triples (rational, integer, integer) encodig the V-filtration index, the weight filtration index and the multiplicity of each spectral pair

!!! note "on planned implementation"
    wrap spectrum und spectral pairs from gmssing.lib, wrap spectrumNND from spectrum.lib

### bernstein

!!! note "Diskussionsbedarf" 
    Das koennen wir heute vermutlich viel besser ueber Viktor's Implementation als mit der gmssing.lib -- mein Vorschlag: wrappen und dokumentieren, aber nicht bei hypersurface singularities, sondern in der Umgebung von Plural"

### vfiltration

!!! note "on planned implementation"
    wrap from gmssing.lib, only include in OSCAR, if Matthias Sch. provides good textbook source for theory
 
### Polar multiplicities and Euler Obstruction

**Wishlist** not available yet
**source** Tajima und Nabeshima

### Le numbers

**Wishlist** not available yet
**source** Cisneros-Le-Seade, volume II

## Classification of Hypersurface Singularities

Arnold's famous classifier, see textbook [AGV85](@cite), volume I, provides an algorithm to determine the type of a given function germ with isolated singularity at the origin w.r.t. R-equivalence up to corank 3, Milnor number 16 and modality 2.  
This classifier and a classifier for hypersurface germs of modality at most 2
and corank at most 2 are available as well as several related classification
tools.  

!!! note "on planned implementation"
    Janko: what do we wrap from where, which one should be named classify_singularity? Do we want one top-level function which then calls the appropriate one via algorithm=...?

### classify

This variant of the classifier provides the type and indices, 
but not the values of the moduli from Arnold's list.  

Provides:** list of string and integer vector

!!! note "on planned implementation"
    wrap by separating pretty printing from returned data in Singular and do the pretty printing in OSCAR

### quickclass

This variant provides a sophisticated guess on the type of the hypersurface singularity via the Mather-Yau theorem, which is often sufficient in smaller practical applications or in combination with additional knowledge on the singularity.  

**Provides:** list of pairs of string and integer vector  

!!! note "on planned implementation"
    wrap by separating pretty printing from returned data in Singular and do the pretty printing in OSCAR

## Deformations and Unfoldings

### Introductory remarks

For a general treatment of deformations of space germs and unfoldings of map germs see the respective sections in [Deformations of Space Germs](@space_germs_deform) and in 
<-- DeformationsOfMapsGerms (HIER ERST DEN LINK SETZEN, WENN MAP GERMS EXISTIERT...) --> A $k$-parameter unfolding of a power series $f \in {\mathbb C}\{x_1,\dots,x_n\}$ is a power series $F \in {\mathbb C}\{x_1,\dots,x_n,t_1,\dots,t_k\}$ such that $F(x,0)=f(x)$, i.e.

```math
F(x,t) = f(x) + \sum_{|a| \geq 1} g_a(x) t^a
```

for suitable $g_a(x) \in {\mathbb C}\{x_1,\dots,x_n\}$.  

Taking the point of view of map germs, an unfolding of $f$ induces an unfolding
$F: ({\mathbb C}^n,0) \times ({\mathbb C}^k,0) \longrightarrow ({\mathbb C},0)$ of the map germ $f: ({\mathbb C}^n,0) \longrightarrow ({\mathbb C},0)$.  

In the hypersurface case, every unfolding of $f$ also induces a deformation of
the space germe $(V(f),0)$ as

```math
(V(f),0) \stackrel{\iota}{\hookrightarrow} (V(F),0)\\ 
\downarrow  \phantom{\;\;\longrightarrow\;\;} \;\; \downarrow \pi\\
\quad\quad \{0\} \;\;\; \longrightarrow\; ({\mathbb C}^k,0)
```

because the flatness condition on $\pi$ holds trivially for hypersurfaces. The
deformation may be understood as a family with base $({\mathbb C}^k,0)$ and special fibre
$(V(f),0)$.

The base of a versal family ( for more information on versality see [Deformations of Space Germs](@space_germs_deform) ) may be chosen to be $({\mathbb C}^\tau,0)$ for a finitely determined hypersurface germ, where $\tau$ is the Tjurina number of the germ. A basis $\{ g_1,\dots ,g_s\}$  of the Tjurina algebra provides the necessary information to write down a versal family as
$(V(F),0) \subset ({\mathbb C}^{n+\tau},0)$ where 

```math
F(x,t_1,\dots,t_{\tau})= f + \sum_{i=1}^{\tau} t_i g_i
``` 

### Tjurina Algebra 
 
The Tjurina algebra of a hypersurface germ can be computed using

```julia 
   tjurina-algebra(f::MPolyElemLoc)
   tjurina-algebra(X::HypersurfaceGerm)
```

**Provides:** MPolyLocQuo

or alternatively (using the fact that the Tjurina algebra is a special 
instance of the $T^1$ of a space germ and that a HypersurfaceGerm is a special instance of a space germ):  

```julia
   T1-module(X::SpaceGerm)
   T1-basis(X::SpaceGerm)
```

!!! note
    Observe that the return value of tjurina_algebra is a localization of a polynomial ring, whereas the one of T1-module is a quotient of a rank 1 free module over a localization of a polynomial ring. For finitely determined hypersurface singularities, both of these objects also carry the structure of a a finite dimensional vector space over the underlying field. A basis thereof is returned by T1-basis.  

!!! note "on planned implementation"
write in OSCAR -- good example for scritping users

### Versal Unfolding

A versal unfolding of a given (computable,) finitely determined  $f \in {\mathcal O}_n$ may be computed using

```julia
   versal_unfolding(f::MPolyLocElem)
```

**Requires:** (V(f),0) finitely determined  

**Provides:** MPolyLocElem in new ring ${\mathcal O}_{n+\tau}$  

!!! note "on planned implementation"
    uses tjurina_algebra

### Versal Deformation

Taking the perspective of space germs, a versal family may be computed
via

```julia
   versal_deformation(X::SpaceGerm)
```

**Provides:** morphisms $\iota$ and $\pi$ of space germs as described above 

!!! note "on planned implementation"
    uses tjurina_algebra

### Kernel of Kodaira-Spencer map/locally trivial subfamilies
