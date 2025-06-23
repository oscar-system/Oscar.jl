# [2025-04-04] FIRST ATTEMPT at defining structs for specifying random params

#####
#####  DON'T WASTE YOUR TIME (& my time) nitpicking about indentation, punctuation, etc. !!!
#####  Also function names can be altered/improved at some later stage...
#####

# The main purpose is to explore possible good public UIs for specifying how
# to generate random values in complicated rings...

# A RandomParams object, RP, offers the following  (similar to Julia Sampler)
# (+) RAND(RP) gives a random value when passed
# (+) RP.counter says how many times a random value has been generated
# (+) RP.num_different_values approximates how many different possible
#     values can be obtained from RP
# (+) extend_range!(RP)  might extend the range of values RP can produce
#     [design still not settled... further changes likely]
# (+) NYI: parent(RP), base_ring(RP), ...

# Naturally questions of semantics arise too!

# Implementation notes:
# (1) It might be worth storing most data members in a Dict{Symbol, Any} so that
#     printing of a RandomParams object is easier to read; of course, having to use
#     type "Any" is a pain...
# (2) Extreme execution efficiency is NOT A CONCERN


abstract type RandomParams end;



function cap_to_Int64(n::T)  where {T <: Nemo.IntegerUnion}
  @req (n >= 0)  "Must be non-negative"
  max = typemax(Int64)
  (n >= max) && return max
  return Int64(n)
end

# used in RAND(::RandomParams_Poly)
function rand_perm(v::Vector{Int})
  n = length(v)
  (n < 2) && return v
  for j in n:-1:2
    k = rand(1:j-1)
    v[j],v[k] = v[k],v[j]  # poor man's swap(v[j],v[k])
  end
  return v
end


# range conversion -- make every range a StepRange{ZZRingElem,ZZRingElem} with positive step
# Julia's UnitRange is too limiting: does not work with ZZRingElem

function range_ZZ(R::ZZRingElemUnitRange)
  return (R.start):ZZ(1):(R.stop)
end

function range_ZZ(R::StepRange{ZZRingElem,ZZRingElem})
  if is_negative(R.step)
    return reverse(R)
  end
  return R;
end

function range_ZZ(R::UnitRange{T})  where {T <: Integer}
  return ZZ(R.start):ZZ(1):ZZ(R.stop)
end

function range_ZZ(R::StepRange{T})  where {T <: Integer}
  if is_negative(R.step)
    return ZZ(R.stop):ZZ(-R.step):ZZ(R.start)
  end
  return ZZ(R.start):ZZ(R.step):ZZ(R.stop)
end


# range doubling -- !!does not double the step size!!

# "semi-smart" doubling of a StepRange
# Looks complicated, but isn't really.
# If range is non-neg, just extend the stop
# If range is non-posg, just extend the start
# O/w straddles zero, so extend both start & stop
# If step is not 1, we keep the residue class fixed (I hope)
function double_range(range::StepRange{ZZRingElem,ZZRingElem})
  # do nothing if range is a singleton
  if range.start == range.stop
    return range
  end
  if range.start == 0
    return (range.start):(range.step):(2*range.stop)
  end
  if range.stop == 0
    return (2*range.start):(range.step):(range.stop)
  end
  width = div(range.stop-range.start, range.step)  # should be exact division
  if range.start > 0
    return (range.start):(range.step):(range.stop + range.step*width)
  end
  if range.stop < 0
    return (range.start - range.step*width):(range.step):(range.stop)
  end
  half_width = div(1+width,2)
  delta = range.step*half_width
  return (range.start-delta):(range.step):(range.stop+delta)
end

function increase_deg_range(range::StepRange{Int,Int})
  return (range.start):(range.step):(range.stop+range.step)  # !! should check for overflow !!
end


# Default: always gives false
function is_compatible(::RandomParams, ::Ring)
  return false
end


function RAND(R::Ring)
  # This should always work if the ring R is finite!
  # Maybe check whether a global flag "PhysicistMode" is set?
  @req (physicist_mode || is_finite(R))  "Please specify a rand_range, or activate physicist mode"
  return RAND(rand_range(R))
end

function rand_range(R::Ring)
  error("Not (yet) defined")
end


#-------------------------------------------------------
# Integers

mutable struct RandomParams_ZZ  <: RandomParams
  range::StepRange{ZZRingElem,ZZRingElem}
  nz::Bool
  counter::Int64  # incremented each time a new random value is generated
  num_different_values::Int64  #

  function RandomParams_ZZ(lwb_upb::StepRange{ZZRingElem,ZZRingElem}, non_zero::Bool)
    @req (lwb_upb.start <= lwb_upb.stop)   "Empty range"
    @req !(non_zero && lwb_upb == 0:0)  "Cannot request non-zero elements from zero range"
    result = new()
    result.range = lwb_upb
    result.nz = non_zero
    result.counter = 0
    result.num_different_values = cap_to_Int64(length(result.range))
    return result
  end
end

function is_compatible(::RandomParams_ZZ, ::ZZRing)
  return true
end

function rand_range(::ZZRing; ZZ_range::AbstractRange{T}=-99:99, non_zero::Bool=false)  where {T <: Nemo.IntegerUnion}
  return RandomParams_ZZ(range_ZZ(ZZ_range), non_zero)
end

# function symmetric_range(::ZZRing, B::ZZRingElem)  # allow more integer types?
#   @req (B >= 0)  "Empty symmetric range"
#   return RandomParams_ZZ(-B:ZZ(1):B, false)
# end

# function nonneg_range(::ZZRing, B::ZZRingElem)  # allow more integer types?
#   @req (B >= 0)  "Empty non-neg range"
#   return RandomParams_ZZ(ZZ(0):ZZ(1):B, false)
# end



function extend_range!(params::RandomParams_ZZ)
  params.range = double_range(params.range)
  ## update counter & num_different_values...
  params.counter = 0 # ? just reset counter ?
  params.num_different_values = cap_to_Int64(length(params.range))
  return params
end



function RAND(params::RandomParams_ZZ)
  params.counter += 1
  val = ZZ(0)  # just for scoping
  # The constructor has checked that this loop cannot be infinite
  while true
    val = ZZ(rand(params.range))
    if (!params.nz || !is_zero(val))
      break;
    end
  end
  return val
end



#-------------------------------------------------------
# Rationals

mutable struct RandomParams_QQ  <: RandomParams
  numer_range::StepRange{ZZRingElem,ZZRingElem} # guaranteed non-empty
  denom_range::StepRange{ZZRingElem,ZZRingElem}  # contains at least 1 non-zero value
  nz::Bool
  counter::Int64  # incremented each time a new random value is generated
  num_different_values::Int64  #

  function RandomParams_QQ(numer::StepRange{ZZRingElem,ZZRingElem}, denom::StepRange{ZZRingElem,ZZRingElem}, non_zero::Bool)
    @req (numer.start <= numer.stop)   "Empty numer range"
    @req (!(denom.start == 0 && denom.stop == 0) && (denom.start <= denom.stop))   "Empty denom range"
    @req !(non_zero && numer == 0:0)  "Cannot request non-zero elements from zero range"
    result = new()
    result.numer_range = numer
    result.denom_range = denom
    result.nz = non_zero
    result.counter = 0
    num_values = div(3*ZZ(length(numer))*ZZ(length(denom)), 5)  # crude approx (asymptotically good estimate)
    result.num_different_values = cap_to_Int64(num_values)
    return result
  end
end

function is_compatible(::RandomParams_QQ, R::QQField)
  return true
end

function rand_range(::QQField;
                    numer_range::AbstractRange{T1}=-99:99,
                    denom_range::AbstractRange{T2}=1:max(1,-numer_range.start,numer_range.stop),
                    non_zero::Bool=false)  where {T1 <: Nemo.IntegerUnion, T2 <: Nemo.IntegerUnion}
  @req (numer_range.start <= numer_range.stop)  "Empty numer range"
  @req !(denom_range.start == 0 && denom_range.stop == 0)  "Empty denom range"
  return RandomParams_QQ(range_ZZ(numer_range), range_ZZ(denom_range), non_zero)
end

# function symmetric_range(::QQField, B::ZZRingElem; non_zero::Bool=false)  # allow more integer types?
#   @req (B >= 0)  "Empty symmetric numer range"
#   return RandomParams_QQ(range_ZZ(-B:B), range_ZZ(1:max(1,B)), non_zero)
# end

# function nonneg_range(::QQField, B::ZZRingElem; non_zero::Bool=false)  # allow more integer types?
#   @req (B > 0)  "Empty non-neg numer range"
#   return RandomParams_QQ(0:B, 1:B, non_zero)
# end



# ??need variants to extend just numer or just denom??
function extend_range!(params::RandomParams_QQ)
  params.numer_range = double_range(params.numer_range)
  params.denom_range = double_range(params.denom_range)
## update counter & num_different_values
  params.counter = 0 # ? just reset counter ?
  num_values = div(3*ZZ(length(params.numer_range))*ZZ(length(params.denom_range)), 5)  # crude approx
  result.num_different_values = cap_to_Int64(num_values)
  return R
end


function RAND(params::RandomParams_QQ)
  params.counter += 1
  # The constructor has checked that this loop cannot be infinite
  numer = ZZ(0) # just for scoping
  denom = ZZ(0) # just for scoping
  while true
    numer = rand(params.numer_range)
    if (!params.nz || !is_zero(numer))
      break;
    end
  end
  while true
    denom = rand(params.denom_range)
    if !is_zero(denom)
      break;
    end
  end
  return QQ(numer//denom)
end


#-------------------------------------------------------
# Finite fields --> TODO

#-------------------------------------------------------
# Finite quotients of ZZ --> TODO


#-------------------------------------------------------
# Univariate polynomials

# QUESTION: what does the deg_range parameter mean?
#    ?ANS?: pick a target deg unif rnd from range, then gen poly of this degree (unless top coeff happens to be 0)

# ??? Allow caller to specify density, num_of_terms (need RandomSubset of 1:n)

mutable struct RandomParams_Poly  <: RandomParams
  poly_ring::PolyRing
  deg_range::StepRange{Int,Int}  # step is positive
  coeff_params::RandomParams
  num_terms::StepRange{Int,Int}  # step is positive  **NOTE**  0:0 means "all"
  nz::Bool
  counter::Int64  # incremented each time a new random value is generated
  num_different_values::Int64  #


  function RandomParams_Poly(P::PolyRing, d::StepRange{Int,Int}, nterms::StepRange{Int,Int}, coeff_range::RandomParams, non_zero::Bool)
    @req (d.start >= 0 && d.stop >= d.start)   "Deg range must be non-negative & non-empty"
    @req ((nterms.start == 0 && nterms.stop == 0) || (nterms.start > 0 && nterms.stop >= nterms.start))  "Number of terms must be positive"
    @req (1+d.start >= nterms.start)  "Impossible combination of deg range and num terms"
    @req is_compatible(coeff_range, coefficient_ring(P))  "Params for coeff ring are not compatible"
    result = new()
    result.poly_ring = P
    result.deg_range = d
    result.coeff_params = coeff_range
    result.num_terms = nterms
    result.nz = non_zero
    result.counter = 0
    max_nterms = (nterms == 0:0) ? length(d) : min(length(d), nterms.stop)
    log_est = max_nterms * log2(coeff_range.num_different_values)
    result.num_different_values = (log_est >= 63) ? typemax(Int64) : floor(Int64, exp2(log_est))
    return result
  end
end

function is_compatible(P::RandomParams_Poly, R::PolyRing)
  return (P.poly_ring == R)
#  return is_compatible(P.coeff_params, coefficient_ring(R))
end

function StepRange_Int(r::StepRange{T1,T2})  where {T1 <: Nemo.IntegerUnion, T2 <: Nemo.IntegerUnion}
  if r.step < 0
    return Int(r.stop):Int(r.step):Int(r.start)
  end
  return Int(r.start):Int(r.step):Int(r.stop)
end

function StepRange_Int(r::UnitRange{T})  where {T <: Nemo.IntegerUnion}
  return Int(r.start):Int(1):Int(r.stop)
end

  

# Way to express prefernce for "extend policy"?
function rand_range(P::PolyRing{T};
                    num_terms::AbstractRange{ZZ1}=0:0 #=0:0 means "all" terms=#,
                    coeff_range::RandomParams=rand_range(coefficient_ring(P)),
                    deg::AbstractRange{ZZ2}=0:4,
                    non_zero=false)  where {T <: RingElem, ZZ1 <: Nemo.IntegerUnion, ZZ2 <: Nemo.IntegerUnion}
  deg = StepRange_Int(deg)
  num_terms = StepRange_Int(num_terms)
  return RandomParams_Poly(P, deg, num_terms, coeff_range, non_zero)
end

function extend_range!(R::RandomParams_Poly)
  R.deg_range = increase_deg_range(R.deg_range)
  R.coeff_params = extend_range!(R.coeff_params);
## update counter & num_different_values
  R.counter = 0 # ? just reset counter ?
  max_nterms = (R.num_terms == 0:0) ? length(R.d) : min(length(R.d), R.num_terms.stop)
  log_est = max_nterms * log2(R.coeff_params.num_different_values)
  result.num_different_values = (log_est >= 63) ? typemax(Int64) : floor(Int64, exp2(log_est))
  return R
end

function RAND(params::RandomParams_Poly)
  params.counter += 1
  x = gen(params.poly_ring)
  target_deg = rand(params.deg_range)
  # use heuristic test for inability to obtain a non-zero result
  f = zero(params.poly_ring)
  for trial in 1:30
    exps = collect(0:target_deg)
    t = rand(params.num_terms)
    if params.num_terms != 0:0 && t < 1+target_deg
      exps = sort(first(rand_perm(exps),t))  # ?sort useful?
    end
    # f is always zero here
    for j in exps
      f += RAND(params.coeff_params)*x^j
    end
    if !(params.nz && is_zero(f))
      break;
    end
  end    
  (params.nz && is_zero(f)) && error("Unable to obtain non-zero result")
  return f
end


#-------------------------------------------------------
# Multivariate polynomials

# Box distrib;  simplex distrib (same as homog?);  subset of a support

# Simplex: see http://blog.geomblog.org/2005/10/sampling-from-simplex.html
#  Intermediate:  L := [ -std::log(rand_U01()) | i in 1..n+1 ];
#  Rescale:       L := (1/sum(L)) * L; // sum is 1, ignore last coord

# Option to make the result homogeneous (wrt weight vector)?


mutable struct RandomParams_MPoly  <: RandomParams
  mpoly_ring::MPolyRing
  deg_ranges::Vector{StepRange{Int,Int}}  # step is positive
  coeff_params::RandomParams
  num_terms::StepRange{Int,Int}  # step is positive  **NOTE**  0:0 means "all"
  nz::Bool
  counter::Int64  # incremented each time a new random value is generated
  num_different_values::Int64  #


  function RandomParams_MPoly(Rxyz::MPolyRing, degs::Vector{StepRange{Int,Int}}, nterms::StepRange{Int,Int}, coeff_range::RandomParams, non_zero::Bool)
    @req (length(degs) == ngens(Rxyz))  "Incompatible length of vector of degree ranges"
    @req all(d -> (d.start >= 0 && d.stop >= d.start), degs)   "Deg ranges must be non-negative & non-empty"
    @req ((nterms.start == 0 && nterms.stop == 0) || (nterms.start > 0 && nterms.stop >= nterms.start))  "Number of terms must be positive"
#???    @req (1+d.start >= nterms.start)  "Impossible combination of deg range and num terms"
    @req is_compatible(coeff_range, coefficient_ring(Rxyz))  "Params for coeff ring are not compatible"
    result = new()
    result.mpoly_ring = Rxyz
    result.deg_ranges = degs
    result.coeff_params = coeff_range
    result.num_terms = nterms
    result.nz = non_zero
    result.counter = 0

    max_nterms = min(length(nterms), cap_to_Int64(prod(length.(degs))))
    log_est = max_nterms * log2(coeff_range.num_different_values)
    result.num_different_values = (log_est >= 63) ? typemax(Int64) : floor(Int64, exp2(log_est))
    return result
  end
end

function is_compatible(P::RandomParams_MPoly, Rxyz::MPolyRing)
  return (P.mpoly_ring == Rxyz)
end



function rand_range(Rxyz::MPolyRing{T};
                    num_terms::AbstractRange{ZZ1}=5:5 #=0:0 means "all" terms=#,
                    coeff_range::RandomParams=rand_range(coefficient_ring(Rxyz)),
                    degs::Vector{RANGE}=[2:3  for _ in 1:ngens(Rxyz)],
                    non_zero=false)  where {T <: RingElem, ZZ1 <: Nemo.IntegerUnion, ZZ2 <: Nemo.IntegerUnion, RANGE <: AbstractRange{ZZ2}}
  degs = StepRange_Int.(degs)
  num_terms = StepRange_Int(num_terms)
  return RandomParams_MPoly(Rxyz, degs, num_terms, coeff_range, non_zero)
end


# THIS IS SURELY WRONG
function extend_range!(R::RandomParams_MPoly)
  R.deg_ranges = (r -> increase_deg_range(r)).(R.deg_ranges)
  R.coeff_params = extend_range!(R.coeff_params);
## update counter & num_different_values
  R.counter = 0  # ? just reset counter ?
  max_nterms = min(length(R.num_terms), cap_to_Int64(prod(length.(R.deg_ranges))))
  log_est = max_nterms * log2(R.coeff_params.num_different_values)
  R.num_different_values = (log_est >= 63) ? typemax(Int64) : floor(Int64, exp2(log_est))
  return R
end


function exponents_to_term(Rxyz::MPolyRing, exps::Vector{Int})
  n = ngens(Rxyz)
  @req  (length(exps) == n)  "Incompatible exponent vector length"
  result = Rxyz(1)
  x = gens(Rxyz)
  for i in 1:n
    result *= x[i]^exps[i]
  end
  return result
end


function RAND(params::RandomParams_MPoly)
  params.counter += 1
  Rxyz = params.mpoly_ring
  x = gens(Rxyz)
  f = Rxyz(0) # just for scoping
  # use heuristic test for inability to obtain a non-zero result
  for trial in 1:30
    # f is always zero here
    for j in 1:rand(params.num_terms)
      f += RAND(params.coeff_params) * exponents_to_term(Rxyz,rand.(params.deg_ranges))  # !!no check for duplicate terms!!
    end
    if !(params.nz && is_zero(f))
      break;
    end
  end    
  (params.nz && is_zero(f)) && error("Unable to obtain non-zero result")
  return f
end

#-------------------------------------------------------
# Matrices (via MatSpace{T})


mutable struct RandomParams_MatSpace  <: RandomParams
  matrix_space::MatSpace
  coeff_params::RandomParams
  density::Float64 # in interval [0,1]
  structure::Symbol  # in [:generic, :symmetric, :anti_symmetric] # is there a better way? # also ?? [:toeplitz, :hankel, :other] ??
  nz::Bool
  counter::Int64  # incremented each time a new random value is generated
  num_different_values::Int64  #


  function RandomParams_MatSpace(S::MatSpace, density::Float64, coeff_range::RandomParams, structure::Symbol, non_zero::Bool)
    @req (is_compatible(coeff_range, base_ring(S)))  "MatSpace and coeff_range are not compatible"
    @req (0 <= density && density <= 1)   "Density must be between 0 and 1"
    @req (structure in [:generic, :symmetric, :anti_symmetric])  "structure name not recognized"
    result = new()
    result.matrix_space = S
    result.coeff_params = coeff_range
    result.density = density
    result.structure = structure
    result.nz = non_zero
    result.counter = 0
    log_coeff_range = log2(coeff_range.num_different_values)
    log_est = nrows(S)*(ncols(S)*log_coeff_range)
    if structure in [:symmetric, :anti_symmetric]
      log_est /= 2
    end
    result.num_different_values = (log_est >= 63) ? typemax(Int64) : floor(Int64, exp2(log_est))
    return result
  end
end

function is_compatible(P::RandomParams_MatSpace, S::MatSpace)
  return (P.matrix_space == S)
end



function rand_range(S::MatSpace{T};
                    coeff_range::RandomParams=rand_range(base_ring(S)),
                    density::Float64=1.0,
                    structure::Symbol=:generic,
                    non_zero=false)  where {T <: RingElem}
  return RandomParams_MatSpace(S, density, coeff_range, structure, non_zero)
end


function extend_range!(R::RandomParams_MatSpace)
  R.coeff_params = extend_range!(R.coeff_params);
## update counter & num_different_values
  R.counter = 0  # ? just reset counter
  log_coeff_range = log2(R.coeff_params.num_different_values)
  log_est = nrows(S)*(ncols(S)*log_coeff_range)
  if R.structure in [:symmetric, :anti_symmetric]
    log_est /= 2
  end
  R.num_different_values = (log_est >= 63) ? typemax(Int64) : floor(Int64, exp2(log_est))
  return R
end



function RAND(params::RandomParams_MatSpace)
  params.counter += 1
  S = params.matrix_space;
  R = base_ring(S)
  M = zero(S)
  structure = params.structure
  # use heuristic test for inability to obtain a non-zero result
  for trial in 1:30
    # M is always zero here
    for i in 1:nrows(M)
      for j in 1:ncols(M)
        if (i == j) && (structure == :anti_symmetric)
          M[i,j] = 0
          continue
        end
        if (i > j) && (structure in [:symmetric, :anti_symmetric])
          if structure == :symmetric
            M[i,j] = M[j,i]
          else
            M[i,j] = -M[j,i]
          end
          continue;
        end
        (rand() < params.density) && continue
        M[i,j] = RAND(params.coeff_params)
      end
    end
    if !(params.nz && is_zero(M))
      break;
    end
  end
  (params.nz && is_zero(M)) && error("Unable to obtain non-zero result")
  return M
end


#=============================================================================
# Example inputs

#-------------------------------------------------------
# Univariate polys

QQx,_ = polynomial_ring(QQ,"x");
RAND(QQx)  # uses default values

r = rand_range(QQx; deg=3:3);
RAND(r)

coeffs = rand_range(QQ; numer_range=-2:2, denom_range=1:1);
r = rand_range(QQx; coeff_range=coeffs);
RAND(r)

# same as above but with degree range specified
r = rand_range(QQx; deg=1:2, coeff_range=coeffs);
L = [RAND(r)  for _ in 1:10000];
count(is_zero, L) # --> 230

# same as above, but with non_zero=true
r = rand_range(QQx; non_zero=true, deg=1:2, coeff_range=coeffs);
L = [RAND(r)  for _ in 1:10000];
count(is_zero, L) # --> 0

# Use a heuristic to avoid infinite loop below:
zero_coeffs=rand_range(QQ; numer_range=0:0, denom_range=1:1);
r = rand_range(QQx; non_zero=true, deg=1:2, coeff_range=zero_coeffs); # impossible request
RAND(r)  # correctly reports it cannot make a non-zero random element

#-------------------------------------------------------
# Multivariate polys

QQxy,_ = polynomial_ring(QQ, ["x", "y"]);

RAND(QQxy)  # default settngs

r = rand_range(QQxy)  # also default settings
RAND(r)

r = rand_range(QQxy; non_zero=true,  coeff_range=rand_range(QQ; numer_range=0:0, denom_range=1:1)) # impossible request
RAND(r)  # correctly reports it cannot make a non-zero random element

# default settings for coeffs...
r = rand_range(QQxy; degs=[1:2, 2:3])
RAND(r)

r = rand_range(QQxy; degs=[1:2, 2:3], num_terms=10:10) # will (almost) surely create duplicate terms
RAND(r)

#-------------------------------------------------------
# Matrix

S = matrix_space(ZZ, 9,9);
r = rand_range(S; density=0.5, structure=:anti_symmetric)
RAND(r)

=============================================================================#
