# Resultants for multivariate polynomial systems

Consider the polynomial system 
```math
F = 
\left\{\begin{array}{ccccc}
f_0 &=& \sum_{\nu \in A_0} a_{0, \nu} \cdot z^\nu &=& 0\\
& \vdots & & & \\
f_m &=& \sum_{\nu \in A_m} a_{m, \nu} \cdot z^\nu &=& 0
\end{array}
\right.
```
over the Laurent polynomial ring ``\mathbb Z[z_1^{\pm 1}, \dots, z_n^{\pm 1}]``.
When ``m \geq n`` the system is in general overdetermined. For the particular 
case ``m = n`` (one more equation than variables) there is a method to determine 
whether or not there is a solution, which applies to most *monomial support sets* 
``\mathbf A = (A_0, \dots, A_n)``. Namely, there exists a hypersurface 
``\nabla \subset \mathbb C^{\mathbf A}`` in the space of coefficients ``\mathbb a``, 
which is the closure of the set of polynomial systems ``F_{\mathbf a}`` for which 
there exists at least one solution in the complex torus ``(\mathbb C^*)^n``:
```math
\nabla = \overline{\left\{\mathbf a \in \mathbb C^{\mathbf A} | \exists \zeta \in (\mathbb C^*)^n : F_{\mathbf a}(\zeta) = 0\right\}}
``` 
The defining equation 
``\Delta \in \mathbb Z[a_{i, \nu} : i = 0,\dots, n, \, \nu \in A_i]`` for ``\nabla`` 
is called the *resultant* for the support set ``\mathbf A``. In general it must be assumed 
to come with a (scheme theoretic) multiplicity. For more details see e.g. [GKZ08](@cite) or [GZ26](@cite). 

We have methods to compute such resultants, going back to ideas of Weyman [Wey94](@cite). 
Namely, for suitable collections of monomial support sets ``\mathbf A = (A_0, \dots, A_n)`` 
as above there exists a toric variety ``X`` without torus factors such that 
```math
\Delta = \det R\pi_*(\mathcal K \otimes \mathcal O(-\alpha))
```
for the projection ``\pi \colon X \times \mathrm{Spec} R \to \mathrm{Spec} R``, 
where ``R = \mathbb Z[a_{i, \nu} : i = 0,\dots, n, \, \nu \in A_i]`` is the ring of 
coefficients and the *tautological polynomials* ``F_i = \sum_{\nu \in A_i} a_{i, \nu} \cdot z^\nu`` give rise to global sections in naturally associated toric line bundles 
on ``X \times \mathrm{Spec} R``. From these one can form a Koszul complex ``K`` 
which sheafifies to a complex of coherent sheaves ``\mathcal K``. The factor 
``\mathcal O(-\alpha)`` is an auxillary twist which can be chosen freely. 

```@docs
    a_resultant(support_sets::Vector{T}; ctx::ResultantCtx, twist) where {T <: Union{Matrix, MatrixElem}}
    a_resultant(F::Vector{T}; ctx::ResultantCtx, twist) where {T<:MPolyRingElem}
```
The above function is mostly for the user's convenience. If one really wants to compute 
resultants, one probably wants to get one's hands on the underlying complex. In fact, 
in many applications the computation of the determinant turns out to be the bottleneck, 
while the complex for the direct image comes out rather quickly. Finding a good twist 
is often key for success.
```@docs
    a_resultant_complex(support_sets::Vector{T}; ctx::ResultantCtx, twist) where {T <: Union{Matrix, MatrixElem}}
    a_resultant_complex(F::Vector{T}; ctx::ResultantCtx, twist) where {T<:MPolyRingElem}
```
A particular application for resultants are *discriminants*. These are varieties in the 
parameter space of families of polynomials, whose fibers are singular hypersurfaces. 
```@docs
    discriminant(f::MPolyRingElem; ctx::ResultantCtx, twist::FinGenAbGroupElem)
    discriminant_complex(f::MPolyRingElem; ctx::ResultantCtx, twist::FinGenAbGroupElem)
    discriminant_complex(A::Union{Matrix, MatrixElem}; tautological_polynomial, ctx::ResultantCtx, twist::FinGenAbGroupElem)
```
