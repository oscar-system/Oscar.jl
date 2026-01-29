```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Schur polynomials

Given a partition $\lambda$ with $n$ parts, the **Schur polynomial** is defined to be
the polynomial

$$s_\lambda := \sum x_1^{m_1}\dots x_n^{m_n}$$

where the sum is taken over all semistandard tableaux $T$ of shape $\lambda$
and $m_i$ is the weight of $i$ in $T$.

There are two different algorithms for the computation of a Schur polynomial
implemented which are automatically selected depending on the size of the
input.

For small integers or if $n\geq 10$, the *combinatorial algorithm* is used.
This algorithm directly applies the above definition.

In the other cases, *Cauchy's bialternant formula*

$$s_\lambda(x_1, \dots, x_n) = \prod_{1\leq i < j \leq n} (x_i - x_j)^{-1}
\begin{vmatrix}
x_1^{\lambda_1 + n - 1} & x_2^{\lambda_1 + n - 1} & \dots & x_n^{\lambda_1 + n - 1} \\
x_1^{\lambda_2 + n - 2} & x_2^{\lambda_2 + n - 2} & \dots & x_n^{\lambda_2 + n - 2} \\
\vdots & \vdots & \ddots & \vdots \\
x_1^{\lambda_n} & x_2^{\lambda_n} & \dots & x_n^{\lambda_n}
\end{vmatrix}$$

is used.
```@docs
schur_polynomial
```
