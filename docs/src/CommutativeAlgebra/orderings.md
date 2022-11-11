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
for the polynomial ring over $C$ in the set of variables $x=\{x_1, \ldots, x_n\}$. Monomials
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

!!! note
    The lexicograpical monomial ordering `lex` specifies the default way of storing and displaying multivariate polynomials in OSCAR (terms are sorted in descending order).
    The other orderings which can be attached to a multivariate polynomial ring are the degree lexicographical ordering `deglex` and the degree reverse lexicographical
	ordering`degrevlex`. Independently of the attached orderings, Gröbner bases can be computed with respect to any monomial ordering. See the section on Gröbner bases.

In this section, we show how to create monomial orderings in OSCAR. After recalling that all monomial orderings can be realized as matrix orderings,
we present a list of orderings which are predefined in OSCAR. Then we discuss weight and block orderings (product orderings). Finally, we address
elimination orderings.

!!! note
    For the convenient construction of block orderings on the set of monomials of a given multivariate polynomial ring, we allow to construct orderings on the
    monomials in blocks of variables, viewing these orderings as partial orderings on the monomials in all variables.

Here are some illustrating examples:

##### Examples

```@repl oscar
S, (w, x) = PolynomialRing(QQ, ["w", "x"])
o = lex([w, x])
canonical_matrix(o)
cmp(o, w^2*x, x^3)
R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
o1 = degrevlex([w, x])
is_global(o1)
canonical_matrix(o1)
o2 = neglex([y, z])
is_local(o2)
canonical_matrix(o2)
o3 = o1*o2
canonical_matrix(o3)
is_mixed(o3)
show(collect(terms((1+w+x+y+z)^2, o3)))
```

## Monomial Comparisons

The `cmp` function should be used for comparing two monomials with regard to a monomial ordering.

```@docs
cmp(ord::MonomialOrdering, a::MPolyElem, b::MPolyElem)
```

Also, the usual iterators for multivariate polynomials have extensions to an
arbitrary ordering.

```@docs
terms(f::MPolyElem, ord::MonomialOrdering)
coefficients(f::MPolyElem, ord::MonomialOrdering)
exponent_vectors(f::MPolyElem, ord::MonomialOrdering)
monomials(f::MPolyElem, ord::MonomialOrdering)
```

## Matrix Orderings

Given a matrix $M\in \text{Mat}(k\times n,\mathbb R)$ of rank $n$, with rows $m_1,\dots,m_k$,
the *matrix ordering* defined by $M$ is obtained by setting
         
$x^\alpha>_M x^\beta  \Leftrightarrow  \;\exists\; 1\leq i\leq k:  m_1\alpha=m_1\beta,\ldots, 
m_{i-1}\alpha\ =m_{i-1}\beta,\ m_i\alpha>m_i\beta$

(here, $\alpha$ and $\beta$ are regarded as column vectors).

!!! note
    By a theorem of Robbiano, every monomial ordering arises as a matrix ordering as above with $M\in \text{GL}(n,\mathbb R)$.

!!! note
    To create matrix orderings, OSCAR allows for matrices with integer coefficients as input matrices.

!!! note
    For orderings such as `lex` and `degrevlex` which are  predefined in OSCAR, using the predefined version is much faster than using a representation as a matrix ordering.

```@docs
matrix_ordering(R::MPolyRing, M::Union{Matrix{T}, MatElem{T}}; check = true) where T
```

As already shown above, OSCAR provides functions to recover defining matrices from given monomial orderings:

```@docs
matrix(ord::MonomialOrdering)
```

```@docs
canonical_matrix(ord::MonomialOrdering)
```

## Predefined Global Orderings

#### The Lexicographical Ordering

The *lexicographical ordering* `lex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i.$

```@docs
lex(R::MPolyRing)
```

#### The Degree Lexicographical Ordering

The *degree lexicographical ordering* `deglex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \deg(x^\alpha) > \deg(x^\beta)  \;\text{ or }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i.$

```@docs
deglex(R::MPolyRing)
```

#### The Reverse Lexicographical Ordering

The *reverse lexicographical ordering* `revlex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i  > \beta_i.$

```@docs
revlex(R::MPolyRing)
```

#### The Degree Reverse Lexicographical Ordering

The *degree reverse lexicographical ordering* `degrevlex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \deg(x^\alpha) > \deg(x^\beta)  \;\text{ or }\;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i.$

```@docs
degrevlex(R::MPolyRing)
```

#### Weighted Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *weighted lexicographical ordering*
`wdeglex(W)`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \text{wdeg}(x^\alpha) > \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i).$

```@docs
wdeglex(R::MPolyRing, W::Vector{Int})
```

#### Weighted Reverse Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *weighted reverse lexicographical ordering*
`wdegrevlex`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \text{wdeg}(x^\alpha) > \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i).$

```@docs
wdegrevlex(R::MPolyRing, W::Vector{Int})
```

## Predefined Local Orderings

#### The Negative Lexicographical Ordering

The *negative lexicographical ordering* `neglex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i < \beta_i.$

```@docs
neglex(R::MPolyRing)
```

#### The Negative Degree Lexicographical Ordering

The *negative degree lexicographical ordering* `negdeglex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \deg(x^\alpha) < \deg(x^\beta)  \;\text{ or }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i.$

```@docs
negdeglex(R::MPolyRing)
```

#### The Negative Reverse Lexicographical Ordering

The *negative reverse lexicographical ordering* `negrevlex` is defined by setting

$x^\alpha > x^\beta \;  \Leftrightarrow \;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i  < \beta_i.$

```@docs
negrevlex(R::MPolyRing)
```

#### The Negative Degree Reverse Lexicographical Ordering

The *negative degree reverse lexicographical ordering* `negdegrevlex` is defined by setting $\;\deg(x^\alpha) = \alpha_1 + \cdots + \alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \deg(x^\alpha) < \deg(x^\beta)  \;\text{ or }\;\exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i.$

```@docs
negdegrevlex(R::MPolyRing)
```

#### Negative Weighted Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *negative weighted lexicographical ordering*
`negwdeglex(W)`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \;  \text{wdeg}(x^\alpha) < \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_1 = \beta_1, \dots, \alpha_{i-1} = \beta_{i-1}, \alpha_i > \beta_i).$

```@docs
negwdeglex(R::MPolyRing, W::Vector{Int})
```

#### Negative Weighted Reverse Lexicographical Orderings

If `W` is a vector of positive integers  $w_1, \dots, w_n$, the corresponding *negative weighted reverse lexicographical ordering*
`negwdegrevlex(W)`  is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha > x^\beta \;  \Leftrightarrow \; \text{wdeg}(x^\alpha) < \text{wdeg}(x^\beta)  \;\text{ or }\;\\
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; \exists \; 1 \leq i \leq n: \alpha_n = \beta_n, \dots, \alpha_{i+1} = \beta_{i+1}, \alpha_i < \beta_i).$

```@docs
negwdegrevlex(R::MPolyRing, W::Vector{Int})
```

## Weight Orderings

If $W$ is a vector of integers  $w_1, \dots, w_n$, and $>$ is a monomial ordering on $\text{Mon}_n(x)$, then
the corresponding *weight ordering* is defined by setting $\;\text{wdeg}(x^\alpha) = w_1\alpha_1 + \cdots + w_n\alpha_n\;$ and

$x^\alpha >_{W} x^\beta \;  \Leftrightarrow \; \text{wdeg}(x^\alpha) > \text{wdeg}(x^\beta)  \;\text{ or }\;
(\text{wdeg}(x^\alpha) = \text{wdeg}(x^\beta)  \;\text{ and }\; x^\alpha > x^\beta).$


```@docs
weight_ordering(W::Vector{Int}, ord::MonomialOrdering)
```

## Block Orderings

The concept of block orderings (product orderings) allows one to construct new monomial orderings from already given ones: If $>_1$ and $>_2$ are monomial orderings on $\text{Mon}_s(x_1, \ldots, x_s)$ and $\text{Mon}_{n-s}(x_{s+1}, \ldots, x_n)$, respectively, then the *block ordering*
$> \; = \; (>_1, >_2)$ on $\text{Mon}_n(x)=\text{Mon}_n(x_1, \ldots, x_n)$ is defined by setting
          
$x^\alpha>x^\beta  \;\Leftrightarrow\;  x_1^{\alpha_1}\cdots x_s^{\alpha_s} >_1 x_1^{\beta_1}\cdots x_s^{\beta_s} \;\text{ or }\;
\bigl(x_1^{\alpha_1}\cdots x_s^{\alpha_s} = x_1^{\beta_1}\cdots x_s^{\beta_s} \text{ and }  x_{s+1}^{\alpha_{s+1}}\cdots x_n^{\alpha_n} >_2
x_{s+1}^{\beta_{s+1}}\cdots x_n^{\beta_n}\bigr).$

!!! note
    The ordering  $(>_1, >_2)$
    - is global (local) iff both $>_1$ and $>_2$ are global (local). Mixed orderings arise by choosing one of $>_1$ and $>_2$ global and the other one local,
    - is an elimination ordering for the first block of variables iff $>_1$ is global.
	
!!! note
    - The definition of a block ordering above subdivides $x$ into a block of initial variables and its complementary block of variables. 
       Block orderings for a subdivision of $x$ into any block of variables and its complementary block are defined similarly and have similar properties.
    - Inductively, one obtains block orderings composed of more than two individual orderings.

In OSCAR, block orderings are obtained by the concatenation of individual  orderings using the `*` operator.

##### Examples

```@repl oscar
R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"])
o = degrevlex([w, x])*degrevlex([y, z])
```

## Elimination Orderings

Let $R = C[x]=C[x_1, \ldots, x_n]$ be a multivariate polynomial ring with coefficient ring $C$.
Fix a subset $\sigma\subset \{1,\dots, n\}$ and write $x_\sigma$  for the set of variables $x_i$ with
$i\in\sigma$. An *elimination ordering for $x\smallsetminus x_\sigma$*  is a monomial ordering
$>$ on $\text{Mon}_n(x)$ which satisfies the following property: If $a$ is a monomial involving one
of the variables in $x\smallsetminus x_\sigma$ , and $b$ is a monomial depending only on the variables in 
$x_\sigma$, then $a > b.$ Computing a Gröbner basis of $I$ with respect to such an ordering provides one way of finding the
intersection $I\cap C[x_\sigma]$, that is, of  *eliminating the variables in $x\smallsetminus x_\sigma$ from $I$*:
The Gröbner basis elements which only depend on the variables in $x_\sigma$ form a Gröbner basis for
$I\cap C[x_\sigma]$ with respect to the restriction of $>$ to the set of monomials in $I\cap C[x_\sigma]$.

!!! note
    The lexicographical ordering is an elimination ordering for each initial set of variables $x_1, \dots, x_k$.
    If only a fixed subset of variables is considered, suitable weight or block orderings as discussed above are more effective.
    The documentation of the `is_elimination_ordering` function below offers examples and non-examples.

## Tests on Monomial Orderings

```@docs
is_elimination_ordering(ord::MonomialOrdering, V::Vector{<:MPolyElem})
```

```@docs
is_global(ord::MonomialOrdering)
```

```@docs
is_local(ord::MonomialOrdering)
```

```@docs
is_mixed(ord::MonomialOrdering)
```


## Module Orderings

Let $R = C[x]=C[x_1, \ldots, x_n]$ be a multivariate polynomial ring with coefficient ring $C$.
Referring to the section on free modules for details, we recall that by a free $R$-module we mean a free
module of type $R^p$ , where we think of $R^p$ as a free module with a given basis, namely the basis
of standard unit vectors. In what follows, $F$ will denote such free $R$-module, and $\{e_1 ,\dots , e_p\}$
will denote the given basis.

A *monomial in $F$*, involving the basis element $e_i$, is a monomial in $R$ times $e_i$.  A term in $F$
is a monomial in $F$ multiplied by a coefficient $c\in C$. Every nonzero element $f\in F$ can 
be uniquely expressed as the sum of finitely many nonzero terms involving  distinct monomials.
These terms (monomials)  are called the *terms (monomials} of $f$.*

A *monomial ordering* on $F$ is a total ordering $>$ on the set of monomials in $F$ such that if
$x^\alpha e_i$ and $x^\beta e_j$ are monomials in $F$, and $x^\gamma$ is a monomial in $R$, then 

$x^\alpha e_i > x^\beta e_j \Longrightarrow x^\gamma x^\alpha e_i > x^\gamma x^\beta e_j.$
In OSCAR, we require in addition that

$x^\alpha e_i > x^\beta e_i \;\text{ iff }\; x^\alpha e_j > x^\beta e_j \;\text{ for all }\; i,j.$

Then $>$ induces a unique monomial ordering on $R$ in the obvious way, and we say that  $>$ is *global*, *local*, or *mixed*
if the induced ordering on $R$ is global, local, or mixed. 

One way of getting a monomial ordering on $F$ is to pick a monomial ordering $>$ on $R$, and extend it to $F$.
For instance, setting

$x^\alpha e_i >  x^\beta e_j \iff x^\alpha > x^\beta \;\text{ or }\; (x^\alpha = x^\beta \;\text{ and }\; i > j)$

gives priority to the monomials in $R$, whereas the ordering defined below gives priority to the components of $F$:

$x^\alpha e_i >  x^\beta e_j \iff i > j \;\text{ or }\; (i = j\;\text{ and } x^\alpha > x^\beta).$

Alternatively, we may wish to use $i < j$ instead of $i > j$ in this definition.

In other words, these orderings are obtained by concatenating a monomial ordering on the monomials of $R$
with a way of ordering the basis vectors of $F$ or vice versa. In OSCAR, we refer to the $i < j$ ordering on the
basis vectors as *lex*, and to the $i > j$ ordering as *revlex*. And, we use the `*` operator for concatenation. 

##### Examples

```@repl oscar
R, (w, x, y, z) = PolynomialRing(QQ, ["w", "x", "y", "z"]);
F = free_module(R, 3)
o1 = degrevlex(R)*revlex(gens(F))
o2 = revlex(gens(F))*degrevlex(R)
```

The induced ordering on the given polynomial ring is recovered as follows:

```@docs
induced_ring_ordering(ord::ModuleOrdering)
```

The comparison function `cmp` as well as the tests `is_global`, `is_local`, and `is_mixed` are also available for module orderings.
