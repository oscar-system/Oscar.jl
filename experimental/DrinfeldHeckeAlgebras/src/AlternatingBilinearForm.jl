#######################################
# Methods for creating and applying alternating bilinear forms
#######################################

#######################################
# Struct and constructors
#######################################

struct AlternatingBilinearForm{T <: RingElem}
  matrix::MatElem{T}

  function AlternatingBilinearForm(m::MatElem{T}) where T <: RingElem
    @req ncols(m) == nrows(m) "Matrix representing an alternating bilinear form must be square"

    @req is_alternating(m) "Input matrix must be skew symmetric with zero diagonal"

    return new{T}(m)
  end
end

alternating_bilinear_form(m::MatElem{T}) where T <: RingElem = AlternatingBilinearForm(m)

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
is_zero(b::AlternatingBilinearForm) = is_zero(matrix(b))
==(a::AlternatingBilinearForm, b::AlternatingBilinearForm) = matrix(a) == matrix(b)
Base.hash(b::AlternatingBilinearForm, h::UInt) = hash(matrix(b), h)
isequal(a::AlternatingBilinearForm, b::AlternatingBilinearForm) = a == b

#######################################
# Application
#######################################

function (b::AlternatingBilinearForm{T})(v::Vector, w::Vector) where T <: RingElem
  B = matrix(b)

  # Check if dimensions match
  n = ncols(B)
  @req length(v) == n && length(w) == n "Arguments must be of dimension $n"

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

  return dot(v_safe, B*w_safe)
end

