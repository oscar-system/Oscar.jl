```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Eigenvalues and -spaces

OSCAR can compute eigenvalues and -spaces exactly and numerically over rings and fields as follows.

## Eigenvalues

Let us define a matrix over the rationals and ask OSCAR to compute its exact eigenvalues.
```jldoctest eigen1
julia> A = QQ[0 -1 0; 1 0 0; 0 0 2]
[0   -1   0]
[1    0   0]
[0    0   2]

julia> eigenvalues(A)
1-element Vector{QQFieldElem}:
 2
```
Note, however, that this only returns the eigenvalues of the matrix over its base field (the rationals in this example). In order to obtain the eigenvalues over a larger field, we can either base change the matrix or pass the desired field as first argument. In the following, we just consider all eigenvalues in the algebraic closure of the rationals.
```jldoctest eigen1
julia> K = algebraic_closure(QQ)
Algebraic closure of rational field

julia> eigenvalues(K, A)
3-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a2: 1.00000*im}
 {a2: -1.00000*im}

julia> A_K = change_base_ring(K, A)
[   {a1: 0}   {a1: -1.00}      {a1: 0}]
[{a1: 1.00}       {a1: 0}      {a1: 0}]
[   {a1: 0}       {a1: 0}   {a1: 2.00}]

julia> eigenvalues(A_K)
3-element Vector{QQBarFieldElem}:
 {a1: 2.00000}
 {a2: 1.00000*im}
 {a2: -1.00000*im}
```
In the following, we ask OSCAR to compute the rounded eigenvalues of a matrix over the complex numbers with a given precision.
```jldoctest
julia> B = AcbField(100)[0 2; 1 0]
[                               0   2.000000000000000000000000000000]
[1.000000000000000000000000000000                                  0]

julia> eigenvalues(B)
2-element Vector{AcbFieldElem}:
 -[1.414213562373095048801688724209 +/- 1.77e-31]
 [1.414213562373095048801688724209 +/- 1.77e-31]
```

```@docs
eigenvalues(M::MatElem{T}) where T <: RingElem
eigenvalues(L::Field, M::MatElem{T}) where T <: RingElem
eigenvalues_simple
eigenvalues_with_multiplicities(M::MatElem{T}) where T <: FieldElem
eigenvalues_with_multiplicities(L::Field, M::MatElem{T}) where T <: RingElem
```

## Eigenspaces

Let us define a matrix over the rationals and ask OSCAR to compute its exact eigenspaces together with its corresponding eigenvalues.
```jldoctest eigen2
julia> A = QQ[0 -1 0; 1 0 0; 0 0 2]
[0   -1   0]
[1    0   0]
[0    0   2]

julia> eigenspaces(A)
Dict{QQFieldElem, QQMatrix} with 1 entry:
  2 => [0 0 1]
```
Note, however, that this only returns the eigenvalues and -spaces of the matrix over its base field (the rationals in this example). In order to obtain the eigenspaces over a larger field, we can either base change the matrix or pass the desired field as first argument. In the following, we just consider all eigenspaces in the algebraic closure of the rationals.
```jldoctest eigen2
julia> K = algebraic_closure(QQ)
Algebraic closure of rational field

julia> eigenspaces(K, A)
Dict{QQBarFieldElem, AbstractAlgebra.Generic.MatSpaceElem{QQBarFieldElem}} with 3 entries:
  {a1: 2.00}     => [{a1: 0} {a1: 0} {a1: 1.00000}]
  {a2: 1.00*im}  => [{a2: -1.00000*im} {a1: 1.00000} {a1: 0}]
  {a2: -1.00*im} => [{a2: 1.00000*im} {a1: 1.00000} {a1: 0}]

julia> A_K = change_base_ring(K, A)
[   {a1: 0}   {a1: -1.00}      {a1: 0}]
[{a1: 1.00}       {a1: 0}      {a1: 0}]
[   {a1: 0}       {a1: 0}   {a1: 2.00}]

julia> eigenspaces(A_K)
Dict{QQBarFieldElem, AbstractAlgebra.Generic.MatSpaceElem{QQBarFieldElem}} with 3 entries:
  {a1: 2.00}     => [{a1: 0} {a1: 0} {a1: 1.00000}]
  {a2: 1.00*im}  => [{a2: -1.00000*im} {a1: 1.00000} {a1: 0}]
  {a2: -1.00*im} => [{a2: 1.00000*im} {a1: 1.00000} {a1: 0}]
  
```



```@docs
eigenspace
eigenspaces
```