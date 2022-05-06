export AbsModuleStdBasis, AbsIdealStdBasis
export odering

abstract type AbsModuleStdBasis{RingElemType} end

abstract type AbsIdealStdBasis{RingElemType} end

### Essential getters #################################################
#
### return the ordering with respect to which the standard 
# basis had been computed
function ordering(G::AbsModuleStdBasis)
end

function ordering(G::AbsIdealStdBasis)
end
  
### return the ambient free module F which contains the elements of G
function ambient_free_module(G::AbsModuleStdBasis)
end

### return the submodule M ⊂ F ≅ Rⁿ of which G is a standard basis, 
# together with its inclusion M ↪ F.
function sub(G::AbsModuleStdBasis; cached::Bool=true)
end

### return the quotient F/M for the submodule M ⊂ F ≅ Rⁿ of which G 
#is a standard basis, together with its projection F →  F/M.
function quo(G::AbsModuleStdBasis; cached::Bool=true)
end

### return a vector of the module elements in the standard basis
function elements(G::AbsModuleStdBasis)
end

### return the leading monomials of `elements(G)`
function leading_monomials(G::AbsModuleStdBasis)
end


### Essential functionality ###########################################

function normal_form(v::FreeModElem, G::AbsModuleStdBasis)
end

### In the case of a finite dimensional quotient F/G with F ≅ Rⁿ 
# a free module over a polynomial ring R = 𝕜[x₁,…,xₙ] with an ordering <, 
# and 𝕜 a field, return a list of elements v₁,…,vᵣ ∈ F which reduce 
# to a 𝕜-basis of F/G.
function k_base(G::AbsModuleStdBasis)
end

### 
# Suppose G is a standard basis for M ⊂ F ≅ Rⁿ for some polynomial ring 
# R = 𝕜[x₁,…,xₙ] over a field 𝕜. Set N = F/M.
#
# Check whether the given weight is compatible with the ordering; if not, complain.
#
# Return a 𝕜-base for the d-th graded slice of N.
function graded_slice(G::AbsModuleStdBasis, d::Int; w::Vector{T}=Vector{Int}()) where {T<:Integer}
end


### Generic functions #################################################

base_ring(G::AbsModuleStdBasis) = base_ring(ambient_free_module(G))

### return the leading module of G
function leading_module(G::AbsModuleStdBasis)
  # can be implemented generically
end

### return the dimension of the support of the module
function dim(G::AbsModuleStdBasis)
end

function is_finite_dimensional(G::AbsModuleStdBasis) 
end

function as_vector_space(G::AbsModuleStdBasis)
  v = k_base(G)
  F = ambient_free_module(G)
  R = base_ring(F)
  kk = coefficient_ring(R) 
  r = length(v)
  V = vector_space(kk, r)
  VectorType = elem_type(V)
  N, p = quo(G)
  ModuleType = elem_type(N)
  function im_func(w::VectorType) 
    return sum([w[i]*v[i] for i in 1:length(w)])
  end
  function pr_func(w::ModuleType)
    w_norm = normal_form(lift(w), G)
    result = zero(V)
    for (c, m) in zip(coefficients(w_norm), monomials(w_norm))
      for i in 1:length(v)
        if m == v[i] 
          result = result + c*V[i]
          break
        end
      end
    end
    return result
  end
  f = MapFromFunc(im_func, pr_func, V, N)
  return V, f
end
