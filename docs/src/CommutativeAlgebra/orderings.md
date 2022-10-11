```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["orderings.md"]
```

# Monomial Orderings

Given a coefficient ring $C$ as in the previous section, we write $C[x]=C[x_1, \ldots, x_n]$
for the polynomial ring over $R$ in the set of variables $x=\{x_1, \ldots, x_n\}$. Monomials
in $x=\{x_1, \ldots, x_n\}$ are written using multi--indices: If $\alpha=(\alpha_1, \ldots, \alpha_n)\in \N^n$,
set $x^\alpha=x_1^{\alpha_1}\cdots x_n^{\alpha_n}$ and

$\text{Mon}_n(x) :=  \text{Mon}(x_1, \ldots, x_n) := \{x^\alpha \mid \alpha \in \N^n\}.$

A *monomial ordering* on $\text{Mon}_n(x)$ is a total  ordering $>$ on $\text{Mon}_n(x)$ such that

$x^\alpha > x^\beta \Longrightarrow x^\gamma x^\alpha > x^\gamma  x^\beta,
\; \text{ for all }\; \alpha, \beta, \gamma \in \mathbb N^n.$

A monomial ordering $>$ on $\text{Mon}_n(x)$ is called
- *global* if $x^\alpha > 1$ for all $\alpha \not = (0, \dots, 0)$,
- *local* if  $x^\alpha < 1$ for all $\alpha \not = (0, \dots, 0)$, and
- *mixed* if it is neither global nor local.

!!! note
    - A monomial ordering on $\text{Mon}_n(x)$ is global iff it is a well-ordering.
    - To give a monomial ordering on $\text{Mon}_n(x)$ means to give a total ordering $>$ on $ \N^n$ such that
	   $\alpha > \beta$ implies $ \gamma + \alpha > \gamma  + \beta$ for all $\alpha , \beta, \gamma \in \N^n.$
       Rather than speaking of a monomial ordering on $\text{Mon}_n(x)$, we may, thus, also speak of a
	   (global, local, mixed) monomial ordering on $\N^n$.
	
Some monomial orderings are predefined in OSCAR. In this section, we fix our notation with respect to  these orderings
and illustrate their creation using several explicit examples. We then discuss block and matrix orderings.

## Predefined Global Orderings

#### The Lexicographical Ordering

The *lexicographical ordering* `lex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i.$

```@docs
lex(V::AbstractVector{<:MPolyElem})
```

#### The Degree Lexicographical Ordering

The *degree lexicographical ordering* `deglex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \deg(x^\alpha) > \deg(x^\beta)  \;\text{ or }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i.$

```@docs
deglex(V::AbstractVector{<:MPolyElem})
```

#### The Reverse Lexicographical Ordering

The *reverse lexicographical ordering* `revlex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i  > \beta_i.$

```@docs
revlex(V::AbstractVector{<:MPolyElem})
```

#### The Degree Reverse Lexicographical Ordering

The *degree reverse lexicographical ordering* `degrevlex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \deg(x^\alpha) > \deg(x^\beta)  \;\text{ or }\;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i.$

```@docs
degrevlex(V::AbstractVector{<:MPolyElem})
```

#### Weighted Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *weighted lexicographical ordering*
`wdeglex(W)`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \text{wdeg}(x^\alpha) > \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i).$

```@docs
wdeglex(V::AbstractVector{<:MPolyElem}, W::Vector{Int})
```

#### Weighted Reverse Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *weighted reverse lexicographical ordering*
`wdegrevlex`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \text{wdeg}(x^\alpha) > \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i).$

```@docs
wdegrevlex(V::AbstractVector{<:MPolyElem}, W::Vector{Int})
```

## Predefined Local Orderings

#### The Negative Lexicographical Ordering

The *negative lexicographical ordering* `neglex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i < \beta_i.$

```@docs
neglex(V::AbstractVector{<:MPolyElem})
```

#### The Negative Degree Lexicographical Ordering

The *negative degree lexicographical ordering* `negdeglex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \deg(x^\alpha) < \deg(x^\beta)  \;\text{ or }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i.$

```@docs
negdeglex(V::AbstractVector{<:MPolyElem})
```

#### The Negative Reverse Lexicographical Ordering

The *negative reverse lexicographical ordering* `negrevlex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i  < \beta_i.$

```@docs
negrevlex(V::AbstractVector{<:MPolyElem})
```

#### The Negative Degree Reverse Lexicographical Ordering

The *negative degree reverse lexicographical ordering* `negdegrevlex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \deg(x^\alpha) < \deg(x^\beta)  \;\text{ or }\;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i.$

```@docs
negdegrevlex(V::AbstractVector{<:MPolyElem})
```

#### Negative Weighted Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *negative weighted lexicographical ordering*
`negwdeglex(W)`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \text{wdeg}(x^\alpha) < \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i).$

```@docs
negwdeglex(V::AbstractVector{<:MPolyElem}, W::Vector{Int})
```

#### Negative Weighted Reverse Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *negative weighted reverse lexicographical ordering*
`negwdegrevlex(W)`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \text{wdeg}(x^\alpha) < \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i).$

```@docs
negwdegrevlex(V::AbstractVector{<:MPolyElem}, W::Vector{Int})
```

## Block Orderings

The concept of block orderings allows one to construct new monomial orderings from already given ones: If $>_1$ and $>_2$ are monomial orderings on $\text{Mon}_s(x_1, \ldots, x_s)$ and $\text{Mon}_{n-s}(x_{s+1}, \ldots, x_n)$, respectively, then the *block ordering*
$>=(>_1, >_2)$ on $\text{Mon}_n(x)=\text{Mon}_n(x_1, \ldots, x_n)$ is defined by setting
          
$x^\alpha>x^\beta  \;\Leftrightarrow\;  x_1^{\alpha_1}\cdots x_s^{\alpha_s} >_1 x_1^{\beta_1}\cdots x_s^{\beta_s} \;\text{ or }\;
\bigl(x_1^{\alpha_1}\cdots x_s^{\alpha_s} = x_1^{\beta_1}\cdots x_s^{\beta_s} \text{ and }  x_{s+1}^{\alpha_{s+1}}\cdots x_n^{\alpha_n} >_2
x_{s+1}^{\beta_{s+1}}\cdots x_n^{\beta_n}\bigr).$
          
Note that $>=(>_1, >_2)$ is global (local) iff both $>_1$ and $>_2$ are global (local). Mixed orderings arise by choosing
one of $>_1$ and $>_2$ global and the other one local.

In Oscar, block orderings are obtained by the concatination of individually given orderings using the `*` operator.

##### Examples

```@repl oscar
R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
o = degrevlex([w, x])*degrevlex([y, z])
```

## Matrix Orderings

Given a matrix $M\in \text{Mat}(k\times n,\mathbb R)$ of rank $n$, with rows $m_1,\dots,m_k$,
the *matrix ordering* defined by $M$ is obtained by setting
         
$x^\alpha>_M x^\beta  \Leftrightarrow  \;\exists\; 1\leq i\leq k:  m_1\alpha=m_1\beta,\ldots, 
m_{i-1}\alpha\ =m_{i-1}\beta,\ m_i\alpha>m_i\beta$

(here, $\alpha$ and $\beta$ are regarded as column vectors).

!!! note
    By a theorem of Robbiano, every monomial ordering arises as a matrix ordering as above with $M\in \text{GL}(n,\mathbb R)$.
    
To create matrix orderings, OSCAR allows for matrices with integer coefficients as input matrices.

```@docs
matrix_ordering(V::AbstractVector{<:MPolyElem}, M::Union{Matrix{T}, MatElem{T}}) where T
```

```@julia
matrix(ord::MonomialOrdering)
```

Return an invertible matrix such that the matrix ordering defined by this matrix is equal to the given ordering `ord`.

##### Example

```@julia
R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
o = degrevlex(R)
matrix(o)
```

## Weight Orderings

If $W$ is a vector of integers  $w_1, \dots, w_n$, and $>$ is a monomial ordering on $\text{Mon}_n(x)$,
the corresponding *weight ordering* is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha >_{W} x^\beta \;  \Leftrightarrow \; \text{wdeg}(x^\alpha) > \text{wdeg}(x^\beta)  \;\text{ or }\;
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; x^\alpha > x^\beta).$


```@julia
weight_ordering(W::Vector{Int}, ord::MonomialOrdering)
```

```@julia
R, (x, y, z) = PolynomialRing(QQ, ["x", "y", "z"]);
W = [1, 0, -1];
o = degrevlex(R)
matrix(o)
oW = weight_ordering(W, o)
matrix(oW)
```

## Tests on Monomial Orderings

```@docs
is_global(ord::MonomialOrdering)
```

```@docs
is_local(ord::MonomialOrdering)
```

```@docs
is_mixed(ord::MonomialOrdering)
```



