```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction
A linear algebraic group is an affine variety with a group structure. These groups can be embedded into a ``GL_n``
and we represent linear algebraic groups as matrix groups. Here we consider reductive linear algebraic groups.

We follow the conventions of [MT11](@cite).

!!! warning
    Currently only simply connected linear algebraic groups of type A_n are supported. In the future, the functionality could be extended to support arbitrary root data.

!!! warning
    Most functionality is currently limited to finite fields.


## Contact

Please direct questions about this part of OSCAR to the following people:
* Janika Peters
* [Max Horn](https://math.rptu.de/en/wgs/agag/people/head/prof-dr-max-horn)

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).
