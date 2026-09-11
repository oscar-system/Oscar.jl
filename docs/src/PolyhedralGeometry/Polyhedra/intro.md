# Introduction

Let $\mathbb{F}$ be an ordered field; the default is that
$\mathbb{F}=\mathbb{Q}$ is the field of rational numbers and other fields are
not yet supported everywhere in the implementation.

A set $P \subseteq \mathbb{F}^n$ is called a *(convex) polyhedron* if it can be
written as the intersection of finitely many closed affine halfspaces in
$\mathbb{F}^n$.  That is, there exists a matrix $A$ and a vector $b$ such that
$$P = P(A,b) = \{ x \in \mathbb{F}^n \mid Ax \leq b\}.$$ Writing $P$ as above
is called an *$H$-representation* of $P$.

When a polyhedron $P \subset \mathbb{F}^n$ is bounded, it is called a *polytope*
and the fundamental theorem of polytopes states that it may be written as the
convex hull of finitely many points.
That is $$P = \textrm{conv}(p_1,\ldots,p_N), p_i \in \mathbb{F}^n.$$
Writing $P$ in this way is called a *$V$-representation*.
Polytopes are necessarily compact, i.e., they form convex bodies.

Each polytope has a unique $V$-representation which is minimal with respect to
inclusion (or cardinality).
Conversely, a polyhedron which is full-dimensional, has a unique minimal
$H$-representation.
If the polyhedron is not full-dimensional, then there is no canonical choice of
an $H$-representation.
