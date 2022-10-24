```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["groebner_bases.md"]
```

# Gröbner Bases

## Computing Gröbner Bases

```@docs
groebner_basis(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool = false)
```
```@docs
std_basis(I::MPolyIdeal, o::MonomialOrdering)
```
See e.g. [GP08](@cite) for the theoretical background on Gröbner- and standard bases.

### Gröbner Bases with transformation matrix

```@docs
groebner_basis_with_transformation_matrix(I::MPolyIdeal; ordering::Symbol = :degrevlex, complete_reduction::Bool=false)
```

    fglm

    Gröbner walks

    Hilbert-driven

!!! warning "Expert functions for Gröbner bases"
    The following functions are low-level implementations of various Gröbner
    basis algorithms with many adjustable arguments. Only use these
    functions directly if you know what you are doing.

```@docs
f4( I::MPolyIdeal; initial_hts::Int=17, nr_thrds::Int=1, max_nr_pairs::Int=0, la_option::Int=2, reduce_gb::Int=1, info_level::Int=0)
```

### Gröbner Bases over the integers

Over the integers the coefficients of the polynomials 
are not invertible, thus their handling when computing
Gröbner bases and normal forms plays an important role. This is done when 
computing strong Gröbner bases which ensure the following property: 
For any element of an ideal its leading term is divisible by a leading term of an 
element of a corresponding strong Gröbner basis.

The textbook [AL94](@cite) provides details on theory and algorithms as well as references.

```@repl oscar
R, (x,y) = PolynomialRing(ZZ, ["x","y"])
I = ideal(R, [2x,3x,4y])
H = groebner_basis(I)
```

## Leading Ideals

```@docs
leading_ideal(g::Vector{T}; ordering::MonomialOrdering) where { T <: MPolyElem }
leading_ideal(I::MPolyIdeal; ordering::MonomialOrdering)
```

## Normal Forms

```@docs
normal_form(f::T, J::MPolyIdeal) where { T <: MPolyElem }
normal_form(A::Vector{T}, J::MPolyIdeal) where { T <: MPolyElem }
```

## Syzygies

```@docs
syzygy_generators(a::Vector{<:MPolyElem})
```

