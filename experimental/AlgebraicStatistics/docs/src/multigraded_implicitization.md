```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Multigraded Implicitization
Let $\phi: S = \mathbb{K}[x_1, \ldots, x_n] \to R = \mathbb{K}[t_1, \ldots, t_m]$ such that $\ker(\phi)$ is homogeneous in a multigrading specified by the columns of a matrix $A \in \mathbb{Z}^{r \times n}$. Then a basis for each homogeneous components of $\ker(\phi)$ can be computed by solving an associated linear system. In particular given a multidegree $\beta = A \alpha$ such that $\alpha \cdot (1, 1, \ldots, 1)  = d$ let $B_\beta$ be a basis for the homogeneous component $S_\beta$. Then a polynomial $f \in S_\beta$ has the form $f = \sum_{m \in B_\beta} c_m m$ and the constraint that $f \in \ker(\phi)$ of course means that $\phi(f) = 0$ which corresponds to a linear system in the variables $c_m$. A basis for the kernel of this linear system corresponds to a basis for $\ker(\phi) \cap S_\beta$. When the multigrading given by $A$ refines total degree this can be done in parallel and a maximal multigrading which can be induced via the map $\phi$ is given by the homogeneity space of the elimination ideal of $\phi$. This algorithm is especially effective when the rank of $ A $ is close to the dimension of $\ker(\phi)$ or when only a subset of generators of low total degree is needed. For more information on this approach to implicitization, see [CH26](@cite). 


```@docs
components_of_kernel(d::Int, phi::MPolyAnyMap; wp::Union{OscarWorkerPool, Nothing}=nothing)
max_grade_domain(phi::MPolyAnyMap)
compute_component(mon_basis::Vector{<:MPolyDecRingElem}, phi::MPolyAnyMap)
jacobian(phi::MPolyAnyMap)
```
