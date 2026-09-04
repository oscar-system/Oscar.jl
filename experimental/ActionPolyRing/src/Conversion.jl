# The content of this file is intended for subroutines of future algorithms

###############################################################################
#
#  Conversion type 
#
###############################################################################

struct __AlgebraicConversionData{
    T <: RingElement,
    A <: ActionPolyRing{T},
    M <: MPolyRing{T},
    F, # Forward map type
    B, # Backward map type
    S  # Supported vars type
}
  apr::A
  mpr::M
  forward_map::F
  backward_map::B
  supported_vars::S
  original_ring_size::Int
end

###############################################################################
#
#  Variable extraction 
#
###############################################################################

function __extract_occ_vars(F::AbstractVector{<:ActionPolyRingElem{T}}) where {T <: RingElement}
  @req !isempty(F) "Cannot extract variables from an empty vector of action polynomials"
  R = parent(first(F))
  @req all(p -> parent(p) === R, F) "All polynomials must belong to the same ring"
  
  var_set = Set{elem_type(R)}()
  for p in F
    union!(var_set, vars(p; sorted=false))
  end
  
  return collect(var_set)
end  

###############################################################################
#
#  Construction 
#
###############################################################################

function __algebraic_conversion_data(apr::ActionPolyRing{T}, jet_vars::AbstractVector{<:ActionPolyRingElem{T}}; revsorted::Bool=true) where {T <: RingElement}
  @req all(is_gen(v) && parent(v) === apr for v in jet_vars) "Invalid jet variables or parent mismatch"
  @req allunique(jet_vars) "The provided vector of jet variables contains duplicates"

  if revsorted
    jet_vars = sort(jet_vars; rev=true)
  end
  
  R = base_ring(base_ring(apr)) 
  n = length(jet_vars)
  S, S_vars = polynomial_ring(coefficient_ring(apr), n)
  
  fwd_images = zeros(S, ngens(R))
  bwd_images = zeros(R, n)
  
  jtu = __jtu_idx(apr) 
  vtj = __vtj(apr) 

  for (i, v) in enumerate(jet_vars)
    raw_idx = jtu[vtj[v]] 
    fwd_images[raw_idx] = S_vars[i]
    bwd_images[i] = gen(R, raw_idx)
  end
  
  fwd_map = hom(R, S, fwd_images)
  bwd_map = hom(S, R, bwd_images)
  
  supported_vars = Set(jet_vars)
  current_size = length(__jtv(apr))
  
  return __AlgebraicConversionData(apr, S, fwd_map, bwd_map, supported_vars, current_size)
end

# A helper for constructing conversion data from a list of polynomials that need not be jet variables. The occurring jet variables
# are sorted automatically, starting with the largest. It is intended to be used for vectors of equations and inequations in algebraic
# systems.
__algebraic_conversion_data_from_polys(apr::ActionPolyRing{T}, F::AbstractVector{<:ActionPolyRingElem{T}}) where {T <: RingElement} =
  __algebraic_conversion_data(apr, __extract_occ_vars(F), revsorted=true)

###############################################################################
#
#  Applying conversions 
#
###############################################################################

# Forward mapping: ActionPolyRingElem -> MPolyRingElem
function (conv::__AlgebraicConversionData{T})(p::ActionPolyRingElem{T}) where {T <: RingElement}
  @req parent(p) === conv.apr "Parent mismatch: The polynomial does not belong to the action polynomial ring of the conversion data"
  @req length(__jtv(conv.apr)) == conv.original_ring_size "The action polynomial ring has generated new jet variables, so the conversion data is outdated"
  @req all(v -> v in conv.supported_vars, vars(p; sorted=false)) "The action polynomial contains variables not present in the conversion data"
  refresh_p = p + zero(parent(p))
  return conv.forward_map(data(data(refresh_p)))
end

# Backward mapping: MPolyRingElem -> ActionPolyRingElem
function (conv::__AlgebraicConversionData{T})(p::MPolyRingElem{T}) where {T <: RingElement}
  @req parent(p) === conv.mpr "Parent mismatch: The polynomial does not belong to the algebraic ring of the conversion data"
  @req length(__jtv(conv.apr)) == conv.original_ring_size "The action polynomial ring has generated new jet variables, so the conversion data is outdated"
  return conv.apr(conv.backward_map(p))
end

(conv::__AlgebraicConversionData{T})(F::AbstractVector{<:ActionPolyRingElem{T}}) where {T <: RingElement} = conv.(F)
(conv::__AlgebraicConversionData{T})(F::AbstractVector{<:MPolyRingElem{T}}) where {T <: RingElement} = conv.(F)

###############################################################################
#
#  IO 
#
###############################################################################

function Base.show(io::IO, ::MIME"text/plain", conv::__AlgebraicConversionData)
  io = pretty(io)
  is_valid = length(__jtv(conv.apr)) == conv.original_ring_size
  
  if !is_valid
    print(io, "Outdated conversion data for ", Lowercase(), conv.apr, "\n")
    print(io, Indent())
    print(io, "Ring has mutated!")
    print(io, Dedent())
  else
    print(io, "Conversion data for ", Lowercase(), conv.apr, " involving variables:\n")
    print(io, Indent())
    original_vars = sort(collect(conv.supported_vars); rev=true)
    join(io, original_vars, ", ")
    print(io, Dedent())
  end
end

function Base.show(io::IO, conv::__AlgebraicConversionData)
  io = pretty(io)
  is_valid = length(__jtv(conv.apr)) == conv.original_ring_size
  
  if is_terse(io)
    print(io, is_valid ? "Conversion data" : "Outdated conversion data")
  else
    print(io, is_valid ? "Conversion data for " : "Outdated conversion data for ")
    print(io, Lowercase(), conv.apr)
  end
end

