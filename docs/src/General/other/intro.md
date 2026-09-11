# Notes for users of other computer algebra systems

This page collects differences that affect users of all such systems.
Notes for users of specific systems follow on separate pages:

- [Notes for GAP users](@ref)
- [Notes for Magma users](@ref)
- [Notes for SageMath users](@ref)
- [Notes for Singular users](@ref)
- [Notes for polymake users](@ref)

!!! note "Help wanted"
    These pages are far from complete.
    If you are missing something here, or if a difference between OSCAR and
    the system you know took you a while to understand,
    please tell us about it, for example by opening an
    [issue on GitHub](https://github.com/oscar-system/Oscar.jl/issues)
    or on [Slack](https://oscar-system.org/slack).
    Contributions to these pages are very welcome.

## General differences

- Julia evaluates `2^100` to `0` because `2` is regarded as a 64 bit integer.
  Write `ZZRingElem(2)^100` to get a long.

- OSCAR makes a subtle but important [distinction between `/` and `//`](@ref subtle_distinction_for_rings).
