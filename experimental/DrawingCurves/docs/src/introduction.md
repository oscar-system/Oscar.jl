```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

This module provides functionality for drawing real plane algebraic curves
defined by a polynomial $f$ in two variables over $\mathbb{Q}$.

```@docs
draw_curve_tikz
```

!!! note "Curves without critical points"
    Currently our code is unable to deal with curves without critical points.
    These are points where both $f$ and its $y$-derivative vanish. Such
    situations arise for example if the curve consists just of a bunch of
    parallel lines, or a parabola, or is empty.
