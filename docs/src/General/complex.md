# Complex Algorithms in OSCAR

On this page we will list some of the more involved algorithmic problems of
OSCAR which you may encounter in the background. For larger examples these may
not terminate, due to lack of memory or time. We will not go into the details
of the respective complexity, as there is sufficient literature.

Often there are several algorithms in OSCAR solving a particular problem, and
trying different alternatives may be worthwile, as some algorithms may not
terminate, while others finish in an instant.

## Groebner and Standard Bases
A standard basis of an ideal is a generating with special properties. A
standard basis is necessary for many (mathematical) low level operations from
commutative algebra.
- Ideal membership
- Radical of an ideal
- Kernel of a ring homomorphism
- Krull dimension of an ideal

## Double Description
A polyhedron may be described as the convex hull of a finite set of points or
as the intersection of finitely many halfspaces. We omit the more complex cases
of unbounded or non-fulldimensional polyhedra here. Computing one description
from the other is done via double description algorithms. Many simple
algorithms on polyhedra need a double description.
- Equality of polyhedra
- Face lattice
- Lattice points
- Hilbert basis
