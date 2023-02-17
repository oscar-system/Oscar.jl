@doc Markdown.doc"""
    markov_basis(S; use_kernel=true)
    markov_basis(P::Polyhedron{fmpq})

# Arguments
- `S::AbstractVector{<:AbstractVector{U}} where U <: Union{Base.Integer, fmpz}`:
Collection of integer vectors that span the columnspace of a matrix M,
 or when `use_kernel=false` collection of integer vectors that span an integer lattice.

- `S::Union{MatElem{U}, <: AbstractMatrix{U}} where U <: Union{Base.Integer, fmpz}`:
A matrix whose integer kernel span an integere Lattice,
or when use_kernel=false a matrix whose rows span an integer lattice.

- `P::Polyhedron{fmpq}`: Polyhedron whose homogenized lattice points span the columnspace of a matrix.

# Returns
A markov basis for the lattice $L$ where $L$ is either the integer kernel of a matrix
when either `use_kernel=true` or the argument is of type `Polyhedron{fmpq}`, otherwise
$L$ is given by the span of the collection of integer vectors `S`.

# Algorithm
Uses the project and lift algorithm see [HM05](@cite) and [DHK12](@cite)

# Examples
```jldoctest
julia> S = [1 0 0 -1; 0 -1 1 0; -1 0 -1 0; 0 1 0 1]
4×4 Matrix{Int64}:
  1   0   0  -1
  0  -1   1   0
 -1   0  -1   0
  0   1   0   1

julia> markov_basis(S)
Warning: unbounded / huge integer variable. Setting  <= (output unimplemented)
1×4 Matrix{fmpz}:
 1  -1  -1  1
```
"""
function markov_basis(S:: Union{MatElem{U}, <: AbstractMatrix{U}};
                      use_kernel=true) where U <: Union{Base.Integer, fmpz}
    S = fmpz_mat(S)
    if !use_kernel
        S = transpose(S)
    end

    pm_M = Polymake.fulton.markov_basis(
        S,
        Polymake.OptionSet(Dict(:use_kernel=>use_kernel))
    )
    return convert(Matrix{fmpz}, pm_M)
end

function markov_basis(S::AbstractVector{<:AbstractVector{U}};
                      use_kernel::Bool=true) where U <: Union{Base.Integer, fmpz}
    T = matrix(ZZ, Vector{Vector{fmpz}}(S))

    if use_kernel
        T = transpose(T)
    end
    return markov_basis(T; use_kernel=use_kernel)
end

function markov_basis(P::Polyhedron{fmpq})
    pm_M = Polymake.fulton.markov_basis(P.pm_polytope)
    return convert(Matrix{fmpz}, pm_M)
end
