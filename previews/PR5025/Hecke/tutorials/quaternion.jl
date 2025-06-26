__HEADER__

using Hecke

# # Quaternion algebras

# ## Creation

Q = quaternion_algebra(QQ, -1, -1)
#----------------------------------------------------------------------------

# Construct the standard basis:

_, i, j, k = basis(Q)
#----------------------------------------------------------------------------

# Verifying the relations:

i^2 == -1 && j^2 == -1 && i * j == k
#----------------------------------------------------------------------------

# Construction of elements:

alpha = 1 + 2*i + 3*j
#----------------------------------------------------------------------------

# Or via directly supplying the coordinates as a vector:

alpha == Q([1, 2, 3, 0])
#----------------------------------------------------------------------------

# This works for also for number fields:

K, sqrt2 = quadratic_field(2)
Q = quaternion_algebra(K, sqrt2, K(3))
alpha = Q([sqrt2, 1, 0, 1])
#----------------------------------------------------------------------------

# ## Properties of elements

# Get the coefficients with respect to the canonical basis:

coefficients(alpha)
#----------------------------------------------------------------------------

# Trace and norm (also reduced version)

tr(alpha), norm(alpha)
#----------------------------------------------------------------------------

trred(alpha), normred(alpha)
#----------------------------------------------------------------------------

# Image of elements under canonical involution:

conjugate(alpha)
#----------------------------------------------------------------------------

normred(alpha) == conjugate(alpha) * alpha
#----------------------------------------------------------------------------

# ## Division

# For division there are the two functions `divexact_left` and `divexact_right`. If `c = divexact_right(a, b)`, then `a == c * b`. So, `divexact_right(a, b)` returns an element `c`, such that `b` becomes a right-divisor of `a`.

_, i, j, k = basis(Q);
#----------------------------------------------------------------------------

divexact_right(k, j)
#----------------------------------------------------------------------------

k == i * j
#----------------------------------------------------------------------------

divexact_left(k, j)
#----------------------------------------------------------------------------

k == j * (-i)
#----------------------------------------------------------------------------

# ## Polynomials

# Polynomials behave very much like polynonomials over commutative rings, except that everything related to divisions needs to specifiy the "side".

Q = quaternion_algebra(QQ, -1, -1)
_, i, j, k = basis(Q)
Qx, x = Q[:x]
f = i * x^2 + j * x
g = i * x
#----------------------------------------------------------------------------

divexact_right(f, g) == x + k
#----------------------------------------------------------------------------

divexact_left(f, g) == x + (- k)
#----------------------------------------------------------------------------

Hecke.divrem_right(f, g)
#----------------------------------------------------------------------------

Hecke.gcd_right(f, g)
#----------------------------------------------------------------------------

# # Splitting of quaternion algebras

Q = quaternion_algebra(QQ, -1, -1)
is_split(Q)
#----------------------------------------------------------------------------

Q = quaternion_algebra(QQ, 1, -1)
is_split(Q)
#----------------------------------------------------------------------------

is_split_with_zero_divisor(Q)
#----------------------------------------------------------------------------

# # Solving norm equations

# Let's solve a norm equation. We want to check whether $2$ is a norm of $\mathbf{Q}(\sqrt{2})$.

K, sqrt2 = quadratic_field(2)
#----------------------------------------------------------------------------

fl, a = is_norm(K, 2);
#----------------------------------------------------------------------------

fl
#----------------------------------------------------------------------------

# Since elements with given norm are in general large, they are represented in special "factored" form:

a
#----------------------------------------------------------------------------

# We can turn this into an ordinary elements using `evaluate`:

b = evaluate(a)
#----------------------------------------------------------------------------

norm(b) == 2
#----------------------------------------------------------------------------

# If we know that a norm equation has a solution, we can directly ask for it:

norm_equation(K, 2)
#----------------------------------------------------------------------------

# # Representation by binary quadratic forms

# Assume that we have two diagonal quadratic forms $q_1 = \langle a_1, a_2 \rangle$ and $q_2 = \langle b_1, b_2 \rangle$ over a field $K$.
# We want to find an element $d$, which is represented both by $q_1$ and $q_2$.

K = QQ;
#----------------------------------------------------------------------------

a1, a2 = 2, 3
#----------------------------------------------------------------------------

b1, b2 = 3, 4
#----------------------------------------------------------------------------

# We form the quadratic form $q = \langle a_1, a_2, -b_1, -b_2\rangle$. Then the task becomes finding an isotropic vector.

q = quadratic_space(K, diagonal_matrix(K, [a1, a2, -b1, b2]));
#----------------------------------------------------------------------------

# Checking whether such an isotropic vector exists:

is_isotropic(q)
#----------------------------------------------------------------------------

fl, v = is_isotropic_with_vector(q)
#----------------------------------------------------------------------------

# To extract the element $d$, we need to evaluate the quadratic form:

d = v[1]^2 * a1 + v[2]^2 * a2
#----------------------------------------------------------------------------

v[1]^2 * a1 + v[2]^2 * a2 == v[3]^2 * b1 + v[4]^2 * b2
#----------------------------------------------------------------------------
