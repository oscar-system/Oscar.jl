#######################################
# Methods for creating and applying alternating bilinear forms
#
# Cassandra Koenen, 2025
#######################################

#######################################
# Struct and constructors
#######################################

struct AlternatingBilinearForm{T <: RingElem}
  matrix::MatElem{T}

  function AlternatingBilinearForm(m::MatElem{T}) where T <: RingElem
    if ncols(m) != nrows(m)
      throw(ArgumentError("Matrix representing an alternating bilinear form must be quadratic")) 
    end
  
    if transpose(m) != -m || !is_zero(diagonal(m))
      throw(ArgumentError("Input matrix must be skew symmetric with zero diagonal")) 
    end
  
    return new{T}(m)
  end
end

const alternating_bilinear_form = AlternatingBilinearForm

#######################################
# Show IO
#######################################

function show(io::IO, b::AlternatingBilinearForm)
  println(io, "Alternating bilinear form, in the standard basis represented by the Gram matrix")
  display(matrix(b))
end

#######################################
# Generic functions
#######################################

matrix(b::AlternatingBilinearForm) = b.matrix
is_zero(b::AlternatingBilinearForm) = is_zero(b.matrix)
==(a::AlternatingBilinearForm, b::AlternatingBilinearForm) = a.matrix == b.matrix
isequal(a::AlternatingBilinearForm, b::AlternatingBilinearForm) = a == b

#######################################
# Application
#######################################

function (b::AlternatingBilinearForm{T})(v::Vector, w::Vector)::T where T <: RingElem
  B = matrix(b)
  
  # Check if dimensions match
  n = ncols(B)
  if length(v) != n || length(w) != n
    throw(ArgumentError("arguments must be of dimension " * string(n))) 
  end
  
  # Check if values are in base ring
  R = base_ring(B)

  v_safe = nothing
  w_safe = nothing
  try
    v_safe = map(x -> R(x), v)
    w_safe = map(x -> R(x), w)
  catch e
    throw(ArgumentError("The given vectors can not be cast into vectors over the base ring."))
  end

  return dot(v, B*w)
end

