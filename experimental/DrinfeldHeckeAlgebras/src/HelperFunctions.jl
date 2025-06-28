################################################################################
# Helper functions for generic Drinfeld-Hecke form generation and validation
#
# Cassandra Koenen, 2025
################################################################################

################################################################################
# Solves the LES Mx = 0 and returns a parametrized solution over the given ring
################################################################################
function solve_and_parametrize(M::MatElem{T}, R::Ring) where {T <: FieldElem}
  nullity, kernel_basis = nullspace(M)
  m = nrows(kernel_basis)
  
  # If nullity = 0, there is only the solution 0
  if nullity == 0
    return fill(R(), m)
  end
  
  # For creating a parametrized solution, we work over the polynomial ring S = R[t1,...tn]
  parameters = nullity == 1 ? ["t"] : ["t" * string(i) for i in 1:nullity]
  S, _ = polynomial_ring(R, parameters)

  # Use kernel basis to create a parametrized solution for a Drinfeld-Hecke form
  sol = fill(S(), m)
  
  for j in 1:nullity
    sol = sol + [S[j] * kernel_basis[i, j] for i in 1:m]
  end

  return sol
end

################################################################################
# Helper function for show methods
################################################################################
function max_column_length(M::MatElem, j::Int)
  result = 0
  
  for i in 1:nrows(M)
    element_length = length(string(M[i,j]))
    
    if element_length > result
      result = element_length
    end
  end

  return result
end
