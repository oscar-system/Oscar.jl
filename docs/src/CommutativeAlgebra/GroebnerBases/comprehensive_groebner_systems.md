# [Comprehensive Gröbner Systems](@id comprehensive_groebner_systems)

Let $K$ be a field, $X=\{x_1,\dots,x_n\}$ a set of variables, and $U=\{u_1,\dots,u_m\}$ a set of parameters.
A comprehensive Gröbner system of an ideal $I\subset K[U][X]$ is a finite list of triples $(E_i, N_i, G_i)$,
with $E_i, N_i\subset K[U]$ and $G_i\subset K[U][X]$, such that

- The union of the segments $V(E_i)\setminus V(N_i)$ covers the parameter space $\overline{K}^m$;
- For every point $a$ of each segment $V(E_i)\setminus V(N_i)$,
  the specialization $G_i(a)$ is a Gröbner basis of the specialized ideal $I(a)\subset\overline{K}[X]$.

In OSCAR, the function below computes a comprehensive Gröbner system.

```@docs
comprehensive_groebner_system
```
