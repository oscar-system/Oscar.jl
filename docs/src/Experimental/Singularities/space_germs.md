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
in a exact smaller field. But unfortunately, also analytic algebras themselves 
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
- [dJP00](@cite)

For the point of view of schemes, we refer to the page on schemes and the references given there; an example of a standard textbook is
- [Har77](@cite).

## Creating Space Germs in OSCAR

In general, space germs in Oscar are created in the following ways:

 * localization of an affine scheme at a point

   ```julia
   SpaceGerm(X::AbsSpec, I::Ideal)
   ```
   where I is a (maximal) ideal describing the chosen ${\mathbb k}$-point on the affine ${\mathbb k}$-scheme X.  
   **Provides:** SpaceGerm.
   ```julia
   germ_at_point(X::AbsSpec, I::Ideal)
   ```
   where I is a (maximal) ideal describing the chosen kk-point on the affine kk-scheme.  

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
also handled in this way in Oscar.

## Basic functionality for space germs

Most of the basic functionality immediately falls back to the underlying 
affine algebra or affine scheme and is just provided on the geometric side 
for convenience and consistence of functionality:

### internal data of a space germ  

 * Pass from the germ $(X,x)$ back to some affine scheme $X$ with the appropriate localization at $x$ by
   ```julia
   representative(X::SpaceGerm)
   ```
   **Provides:** AffineScheme   

 * Given the space germ $(X,x)$, the point $x$ is returned by:

   ```julia
   point(X:SpaceGerm)
   ```
   **Provides:** Ideal

!!! note
    The returned ideal is a prime ideal in the ring of a representative of the germ. At first use of it or at the latest upon the first call of representative, the respective affine scheme is cached and subsequently used for all further purposes requiring a representative

 
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
   is_canonically_isomorphic(X::SpaceGerm, Y::SpaceGerm)
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
     ```julia
     intersect(X::SpaceGerm,Y::SpaceGerm)
     ```
     Computes the intersection of two subgerms of a common larger germ.

     **Provides**: SpaceGerm   

### singular locus

```julia
singular_locus(X::SpaceGerm)
```

   Computes the germ of the singular locus of the given germ.

   **Provides:** SpaceGerm  
  
If a test for smoothness is the goal, the following functions can be used

```julia
is_smooth(X::SpaceGerm)
```

   **Provides:** Boolean Value

<!-- ```julia is_regular(X::SpaceGerm)```  still needs to be implemented -->


<!-- ...to be continued -->
