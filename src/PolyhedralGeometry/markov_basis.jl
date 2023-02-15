function markov_basis(S::Union{<:AbstractVector{<:AbstractVector{U}},
                               MatElem{U},
                               <: AbstractMatrix{U}}, use_kernel=true) where U <: Union{Base.Integer, fmpz}
    hom_S = homogenized_matrix(S, 1)
    return Polymake.fulton.markov_basis(Polymake.Matrix{Polymake.Integer}(hom_S), Polymake.OptionSet(Dict(:use_kernel=>use_kernel)))
end
