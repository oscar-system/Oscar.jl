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

Hypersurface singularities are space germs of the form $(V(f),p)$ 
for some power series f at a point $p$. As explaind in
section [Generalities on Space Germs](@ref space_germ_generalities), the need for exact computations forces OSCAR
to restrict all considerations to the subring of those representable by means 
of an element of the localization of an affine scheme at a point. Where no other reference is specified, we refer the reader to the textbook [GLS07](@cite).

For simplicity of notation, $O_{n,p}$ will denote the localization of the
polynomial ring in $n$ variables at the point $p$ here and $m_p$ denotes its 
maximal ideal. If $p$ is origin, we often omit $p$ in the index.

**Right** **equivalence** of hypersurface singularities takes the point of 
view of the a map f: (CC^n,p) to (CC,0) defined by the power series $f(x)$ 
and allows coordinate changes in the source,i.e.
```math
f \sim_R g :\Longleftrightarrow \exists \varphi \in Aut(O_{n,p}): f = g \circ \varphi.
```
**Contact** **equivalence** takes the point of view of the vanishing locus of f and
allows isomorphism of space germs, i.e. 
```math
f \sim_C g :\Longleftrightarrow (V(f),x) \cong (V(g),x).
```
 
Due to these two perspectives, the majority of functionality for hypersurfaces 
is available both for elements of $O_{n,0}$ and for space germs at an arbitrary
point $p \in CC^n$ with coordinates in a computable subfields of CC. To keep the 
notation simple and cover both cases simultaneously, the following
discussion will be focusing on germs at the origin. However, the functions 
accepting a space germ as input may be called with a hypersurface germ at any  point $p$ subject to the restriction in the previous sentence.

## Basic properties of Hypersurface Singularities

### isolatedness and singular locus

A hypersurface singularity (V(f),0) is called isolated, if it posseses a 
representative on a sufficiently small neighbourhood of 0, whose singular
locus is the origin. By abuse of notation, we also refer to the germ of
this singular locus as the singular locus of the original germ.

```julia
   is_isolated(f::MPolyElemLoc)
   is_isolated(X::SpaceGerm)
```

**Provides:** Boolean value

```julia
   singular_locus(f::MPolyElemLoc)
   singular_locus(X::SpaceGerm)
```

**Provides:** MPolyIdealLoc ???

!!! note "Bemerkung zur Implementation: wrappe bzw. nutze slocus bzw. dim(slocus)"

### Finite determinacy

A hypersurface germ $(V(f),0)$ is called finitely R-determined (or finitely
C-determined), if 
```math
f \equiv g \;\; mod \;\; m^{k+1} \Longrightarrow f \sim_R g \;\;(or \;\; f \sim_C g \;\;respectively).
```

```julia
is_finitely-determined(f::MPolyElemLoc)
is_finitely-determined(X::SpaceGerm)
```

BOOL -- bound separat...

If the germ is finitely determined, the return value gives a (in general non-optimal) determinacy bound; otherwise, it returns the negative integer -1.

**Provides:** integer

!!! note "wir verwenden $\mu+1$, auch wenn sich bessere Schranken ueber $m^{k+2} \subset m^2 Jf$ oder über die höchste Ecke berechnen lassen -- soll der Nutzer die Moeglichkeit haben, eine bessere Schranke konkret anzufordern?"

### Quasihomogeneity

(V(f),0) is called quasihomogeneous, if there are an integer $d$ and a tuple
of positive rational numbers $w = (w_1,\dots,w_n)$ such that all monomials in f
are of $w$-weighted degree $d$, i.e.
$\sum_{i=1}^n w_ia_i =d$  for any monomial $x^a=x_1^{a_1}\cdots x_n^{a^n}$ in the 
support of f. We will refer to these data as a weight system $(d,w)$.

```julia
is_quasihomogeneous(f::MPolyElemLoc)
is_quasihomogeneous(X::SpaceGerm)
```

BOOL -- bound separat...

If the germ is quasihomogeneous w.r.t. a weight system $(d,w)$, this tuple is
returned. Otherwise, a tuple $(0,0,\dots,0)$ is returned.

**Provides:** List of n+1 rational numbers

### is_Morse-critical-point

A germ $(V(f),0)$ is called a Morse critical point, if it is singular and
the Hessian Matrix $H(f)$ of $f$ at $0$ has full rank.

```julia
   is_Morse-critical-point(f::MPolyElemLoc)
   is_Morse-critical-point(X::SpaceGerm)
```

**Provides:** Boolean Value

## Newton polytope and related data

### Newton polytope

For a multivariate power series $f = \sum_{a \in NN^n}$ c_a x^a$, the Newton
polytope is defined as

```math
\Delta (f) = Conv \{a \in NN^n \mid c_a \neq 0\}
```

!!! note  
    This command is only provided for f and not for space germs as the construction depends on the chosen system of parameters.


```julia
   Newton_polytope(f::MPolyElemLoc)
```

**Provides:** ???

!!! note "**Diskussionsbedarf**"
    wie soll die Ausgabe gestaltet sein?  
      A) Polymake Polytope  
      B) Liste der von Polymake errechneten Ecken des Polytopes als natives OSCAR Objekt?    
    (Ich bin fuer ersteres, da wir mit dem Uebergang zu Polytopen in die Polymake-Welt eingetaucht sind)                   

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

**Provides:** ???

!!! note "**Diskussionsbedarf**: wie soll die Ausgabe gestaltet sein?"  
    A) Polymake Polytope  
    B) Liste der von Polymake errechneten Ecken des Polytopes als natives OSCAR Objekt?  
    C) Polymake Polytope plus Zerlegung in Facetten  
    D) Liste der Facetten als MPolyElemLoc mit Nenner 1 (es sind qh-Polynome)  
    (Ich bin fuer D, da man in der Regel sehr schnell auf die einzelnen Facetten einschraenken moechte und diese wieder in OSCAR/Singular-verdaulicher Form braucht)     

### is_convenient

A power series f is called convenient, if its Newton Diagram meets all the
coordinate axes. The implemented algorithm does not explicitly compute the 
Newton diagram!

```julia
   is_convenient(f::MPolyElemLoc)
```

**Provides:** Boolean value

!!! note "**Implementationshinweis** "
    Das Newton Diagram ist hier overkill, das sieht man direkt an der Potenzreihe bzw. an Zähler und Nenner des MPolyElemLoc.

### is_Newton-non-degenerate

Let $f = \sum_{a \in NN^n}$ c_a x^a$ be a multivariate, convenient power series and let $\Gamma(f)$ be its Newton-diagram. For each face $\sigma \in \Gamma(f)$ we define a quasihomogeneous polynomial 
```math
   f_{\sigma}=\sum_{a \in \sigma} c_a x^a
```
f is called Newton-non-degenerate, if for all faces
$f_{\sigma} \in \Gamma(f)$ the germ $(V(f_\sigma),0)$ is non-singular outside 
of $(V(\prod_{i=1}^n x_i),0)$.

```julia
   is_Newton-non-degenerate(f::MpolyElemLoc)
```

**Provides:** Boolean Value

## Invariants of Hypersurface Singularities

Several geometrical, topological and deformation theoretical properties of hypersurface singularities can be described by numerical invariants. Those which are accessible to computation with OSCAR are listed here

### corank

For a singular space germ $(V(f),0) \subset ({\mathbb C}^n,0)$, the corank is the integer
$n-rank(H(f))$ where $H(f)$ denotes the Hessian matrix of $f$, i.e. the $n \times n$ matrix of second order partial derivatives.  

By the Generalized Morse Lemma the following holds for a space
germ $(V(f),0) \subset ({\mathbb C}^n,0)$ of corank $c$
```math
   f \sim_R g(x_1,\dots,x_c) + \sum_{i=c+1}^n x_n^2
```
where $g \in \langle x_1,\dots,x_c \rangle^3$.  

```julia
   corank(f::MPolyElemLoc)
   corank(X::SpaceGerm)
```

**Provides:** integer

!!! note "Dreizeiler, also Nachprogrammieren!"

### multiplicity

For a hypersurface germ $(V(f),0)$ the multiplicity at the origin is the
order of the power series expansion of $f$, i.e.
```math
   mult(V(f),0) = ord_0(f) = min\{k\in {\mathbb N} \mid f \not\in m^{k+1}\}
```

```julia
   multiplicity(f::MPolyElemLoc)
   multiplicity(X::MPolyElemLoc)
``` 

**Provides:** integer

### Milnor number and Milnor Algebra

For an isolated singularity $(V(f),0) \subset (CC,0)$, its Milnor 
fibre has the topological type of a bouquet of $\mu$ $(n-1)$-spheres. This
number $\mu$ is called the Milnor number of $(V(f),0)$ and can be computed
as the vector space dimension of the Milnor algebra:  

```math
\mu(V(f),0) = dim_{CC} (O_n / \langle \frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

For an affine hypersurface $V(f)$ possessing at most isolated singularities, 
the sum over the Milnor numbers at all singular points can be computed as  

```math
\mu(V(f)) = dim_{CC} (CC[x_1,\dots,x_n]  / 
             \langle \frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

```julia
   milnor_number(f::MPolyElemLoc)
   milnor_number(X::SpaceGerm)
   milnor_number(f::MPolyElem)
   milnor_number(X::AffineScheme)
```

**Provides:** integer


```julia
   milnor_algebra(f::MPolyElemLoc)
   milnor_algebra(X::SpaceGerm)
   milnor_algebra(f::MPolyElem)
   milnor_algebra(X::AffineScheme)
``` 

**Provides:** MPolyQuo or MPolyLocQuo (depending on input)

!!! warning  
    In contrast to other functions for space germs, milnor_number and 
    milnor_algebra will not throw an error when called with an affine 
    scheme or a polynomal ring element, but compute a different value!

!!! note "in OSCAR schreiben, gute Tests fuer Uebergang zu endl. Dim VR" 

### Tjurina number

For an isolated singularity $(V(f),0) \subset ({CC}^n,0)$, the 
base of a miniversal family has dimension $\tau$. This number $\tau$ is 
called the Tjurina number of $(V(f),0)$ and can be computed
as the vector space dimension of the Tjurina algebra:

```math
   \tau(V(f),0) = dim_{CC} (O_n / \langle f,\frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

For an affine hypersurface $V(f)$ possessing at most isolated singularities, 
the sum over the Tjurina numbers at all singular points can be computed as

```math
   \tau(V(f)) = dim_{CC} ({CC}[x_1,\dots,x_n]  / 
             \langle f,\frac{\partial f}{\partial x_1},\dots,
                        \frac{\partial f}{\partial x_n}) \rangle
```

```julia
   tjurina_number(f::MPolyElemLoc)
   tjurina_number(X::SpaceGerm)
   tjurina_number(f::MPolyElem)
   tjurina_number(X::AffineScheme)
``` 

**Provides:** integer

!!! warning  
    In contrast to other functions for space germs, tjurina_number will not throw an error when called with an affine scheme or a polynomal ring element, but compute a different value.

!!! note "in OSCAR schreiben, Dreizeiler" 
 
### monodromy

For a detailed introduction to the monodromy of a hypersurface singularity we
refer the reader to [CLS20](@cite), chapter8. Here we assume familiarity 
with the notion and explain what data is provided by the respective function
in OSCAR.  

Denote by $X_t$ a Milnor fibre of a $n$-dimensional isolated hypersurface 
singularity $(X,0) \subset (CC^{n+1},0)$. Then the Monodromy Theorem states
that the eigenvalues of the monodromy $M \in Aut(H^n(X_t),ZZ)$ are roots of
unity and its Jordan blocks are of size at most $(n+1) \timex (x+1)$. The
eigenvalues and the sizes and multiplicities of the corresponding Jordan 
blocks can be determined using OSCAR:

```julia
   monodromy(f::MPolyElemLoc)
   monodromy(X::SpaceGerm)
```

**Provides:** ???

!!! note "Wie stellen wir den Rueckgabewert vernuenftig dar? Bisher sind es drei Listen, besser waere eine Liste von Tupeln mit Eigenwert, Blockgroesse, Vielfachheit, oder nicht?"

!!! note  
    The eigenvalues are not returned as roots of unity, but as a rational
    number $q \in QQ$ such that $e^{2 \pi i q}$ is the desired eigenvalue. 

!!! note "aus gmssing.lib uebernehmen, mondromy.lib ist aelter, weniger allgemein und in bestimmten Faellen instabiler"

### spectrum and spectral pairs

For a detailed treatment of the spectrum of a singularity and of the spectral
pairs, we refer the reader to AGV Band 2 HIER ZITIEREN SOBALD EINTRAG GEMACHT.
It is a powerful invariant, which permits to exclude equivalence or adjacency
of given hypersurface germs in many practical cases.  

In OSCAR the singularita spectrum can be computed using two different methods:


```julia
   spectrum(f::MPolyElemLoc)
   spectrum(X::SpaceGerm)
   spectrumNND(f:MPolyElemLoc)
```

 **Provides:** list of pairs (rational, integer) encoding the spectral numbers and their multiplicities  

!!! note  
    spectrumNND is only applicable, if the input $f$ is Newton non-degenerate. This is checked at the start of the procedure and an error is returned, if $f$ does not have this property. If spectrumNND is applicable, it uses a combinatorial algorithm and is often faster than the function spectrum which relies on the Brieskorn lattice/the Gauss-Manin connection.

```julia
   spectral-pairs(f::MpolyElemLoc)
   spectral-pairs(X::SpaceGerm)
```

**Provides:** list of triples (rational, integer, integer) encodig the V-filtration index, the weight filtration index and the multiplicity of each spectral pair

!!! note "spectrum und spectral pairs aus gmssing.lib, spectrumNND aus spectrum.lib"

### bernstein

!!! note "Diskussionsbedarf: Das koennen wir heute viel besser ueber Viktor's Implementation als mit der gmssing.lib -- mein Vorschlag: wrappen und dokumentieren, aber nicht bei hypersurface singularities, sondern in der Umgebung von Plural"

### vfiltration

!!! note "Anne: mit Matthias Sch. nochmal ueber gute Quelle reden... ohne Quelle: raus"
 
### Polar multiplicities and Euler Obstruction

Nicht vorhanden in Singular!!! sollte man als Arbeit ausgeben
Quelle: Artikel von Tajima und Nabeshima


### Le numbers

Nicht vorhanden in Singular!!! sollte man als Arbeit ausgeben
Quelle: Kapitel in Cisneros-Le-Seade Band II

## Classification of Hypersurface Singularities

Arnold's famous classifier ...AVG Band 1 zitieren... provides an algorithm to determine
the type of a given function germ with isolated singularity at the origin
w.r.t. R-equivalence up to corank 3, Milnor number 16 and modality 2.  
This classifier and a classifier for hypersurface germs of modality at most 2
and corank at most 2 are available as well as several related classification
tools.  

!!! note "Janko, hier hast Du die aktuellere Uebersicht. Sag an, was reingehoert und ich schreibe die Hilfeeintraege entsprechend. "

### classify

This variant of the classifier provides the type as a string and the indices, 
but not the values of the moduli from Arnold's list.  
 
Vorgehen: wrappen

Provides:** list of string and integer vector

!!! note "produziert pretty printing -- wie geht man damit um?"

### quickclass

This variant provides a sophisticated guess on the type of the hypersurface singularity via the Mather-Yau theorem, which is often sufficient in smaller practical applications.  

Vorgehen: wrappen  

**Provides:** list of pairs of string and integer vector  

Bemerkung: produziert pretty printing -- sollte weiterhin optional möglich sein

!!! note "produziert pretty printing -- wie geht man damit um?"

### welche weiteren??

## Deformations and Unfoldings

### Introductory remarks

For a general treatment of deformations of space germs and unfoldings of map germs see the respective sections in [Deformations of Space Germs](@space_germs_deform) and in DeformationsOfMapsGerms (HIER ERST DEN LINK SETZEN, WENN MAP GERMS EXISTIERT...) A $k$-parameter unfolding of a $f \in CC\{x_1,\dots,x_n\}$ is a power series $F \in CC\{x_1,\dots,x_n,t_1,\dots,t_k\}$ such that $F(x,0)=f(x)$, i.e.

```math
F(x,t) = f(x) + \sum_{|a| \geq 1} g_a(x) t^a
```

for suitable $g_a(x) \in CC\{x_1,\dots,x_n\}$.  

Taking the point of view of map germs, an unfolding of $f$ induces an unfolding
$F: (CC^n,0) \times (CC^k,0) \longrightarrow (CC,0)$ of the map germ $f: (CC^n,0) \longrightarrow (CC,0)$.  

In the hypersurface case, every unfolding of $f$ also induces a deformation of
the space germe $(V(f),0)$ as

```math
(V(f),0) \stackrel{\iota}{\hookrightarrow} (V(F),0)\\ 
\downarrow  \phantom{\;\;\longrightarrow\;\;} \;\; \downarrow \pi\\
\quad\quad \{0\} \;\;\; \longrightarrow\; (CC^k,0)
```

because the flatness condition on $\pi$ holds trivially for hypersurfaces. The
deformation my be understood as a family with base $(CC^k,0)$ and special fibre
$(V(f),0)$.

The base of a versal family ( for more information on versality see [Deformations of Space Germs](@space_germs_deform) ) may be chosen to be 
$(CC^\tau,0)$ for a finitely determined hypersurface germ, where $\tau$ is 
the Tjurina number of the germ, i.e. the vector space dimension of its 
Tjurina algebra. A basis $\{ g_1,\dots ,g_s\}$  of the Tjurina algebra 
provides the necessary information to write down a versal family as
$(V(F),0) \subset (CC^{n+\tau},0)$ where 

```math
F(x,t_1,\dots,t_{\tau})= f + \sum_{i=1}^{\tau} t_i g_i
``` 

### Tjurina Algebra
 
The Tjurina algebra of a hypersurface germ can be computed using

```julia 
   tjurina-algebra(f::MPolyElemLoc)
   tjurina-algebra(X::SpaceGerm)
```

**Provides:** MPolyLocQuo

or alternatively (using the fact that the Tjurina algebra is a special 
instance of the T^1 of a space germ):  

```julia
   T1-module(X::SpaceGerm)
   T1-basis(X::SpaceGerm)
```

!!! note
    Observe that the return value of tjurina_algebra is a localization of a polynomial ring, whereas the one of T1-module is a quotient of a rank 1 free module over a localization of a polynomial ring. For finitely determined hypersurface singularities, both of these objects also carry the structure of a a finite dimensional vector space over the underlying field. A basis thereof is returned by T1-basis.  

!!! note "tjurina_algebra ist kurz und liefert ein gutes Beispiel auch fuer proggrammierende Nutzer, die mal spicken wollen -- daher in OSCAR schreiben; T1 siehe Kommentare in spacegerms.md"

### Versal Unfolding

A versal unfolding of a given (computable,) finitely determined  $f \in O_n$ may be computed using

```julia
   versal_unfolding(f::MPolyLocElem)
```

**Requires:** (V(f),0) finitely determined  

**Provides:** MPolyLocElem in new ring $O_{n+\tau}$  

!!! note "In jedem Fall auf der Basis von tjurina_algebra neu schreiben (dreizeiler). Wie soll hier die Wahl der Parameternamen stattfinden? Soll der User den Wortstamm der neuen Variablen angeben?"

### Versal Deformation

Taking the perspective of space germs, a versal family may be computed
via

```julia
   versal_deformation(X::SpaceGerm)
```

**Provides:** morphisms $\iota$ and $\pi$ of space germs as described above at the beginning of the section

!!! note "Kommentare zu Implementation in spacegerms.md beachten. Diskussionsbedarf: Ist das so vernuenftig? Mathematisch ist das genau der Datensatz einer Deformation, aber will man so damit arbeiten?"
