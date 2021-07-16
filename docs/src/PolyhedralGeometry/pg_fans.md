```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["pg_fans.md"]
```

# Polyhedral Fans

## Introduction

Let $\mathbb{F}$ be an ordered field; the most prominent case here is where $\mathbb{F}=\mathbb{Q}$ are the rational numbers.

A nonempty finite collection $\mathcal{F}$ of (polyhedral) cones in $\mathbb{F}^n$, for $n$ fixed, is a (polyhedral) fan if

- the set $\mathcal{F}$ is closed with respect to taking faces and
- if $C,D\in\mathcal{F}$ then $C\cap D$ is a face of both, $C$ and $D$.

## Construction
