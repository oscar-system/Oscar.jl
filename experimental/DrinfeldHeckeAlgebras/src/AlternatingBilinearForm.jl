################################################################################
# Alternating bilinear form
################################################################################

function alternating_bilinear_form(m::MatElem)
  if transpose(m) != -m
    throw(ArgumentError("input matrix must be skew symmetric")) 
  end

  return BilinearForm(m)
end

################################################################################
# Bilinear form
################################################################################

struct BilinearForm{T <: RingElem}
  matrix::MatElem{T}

  function BilinearForm(m::MatElem{T}) where T <: RingElem
    if ncols(m) != nrows(m)
      throw(ArgumentError("matrix defining image of bilinear form must be quadratic")) 
    end
  
    return new{T}(m)
  end
end

function show(io::IO, b::BilinearForm)
  print(io, "Bilinear form defined by image matrix ")
  show(io, b.matrix)
end

################################################################################
# Generic functions
################################################################################

matrix(b::BilinearForm) = b.matrix
is_zero(b::BilinearForm) = is_zero(b.matrix)
==(a::BilinearForm, b::BilinearForm) = a.matrix == b.matrix
isequal(a::BilinearForm, b::BilinearForm) = a == b

################################################################################
# Application
################################################################################

function (b::BilinearForm{T})(v::Vector{T}, w::Vector{T}) where T <: RingElem
  B = b.matrix
  n = ncols(B)
  
  if length(v) != n || length(w) != n
    throw(ArgumentError("arguments must be of dimension " * string(n))) 
  end

  return dot(v, B*w)
end

