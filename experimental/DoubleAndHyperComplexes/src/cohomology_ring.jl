@doc raw"""
    DGAlgCohRing{T} <: NCRing

A struct for cohomology rings `A` which arise from a cochain complex `C` 
which is equipped with a structure of a differential graded (DG) algebra 
by implementing a method for the internal function `mul_cochains`.
"""
mutable struct DGAlgCohRing{T} <: NCRing
  C::SimplicialCochainComplex
  graded_parts::Vector{SubquoModule{T}}
  volume_form::Tuple{Int, SubquoModuleElem{T}}
  vol_form_inc::SubQuoHom
  small_gens::Dict{Int, Vector{SubquoModuleElem{T}}}

@doc raw"""
    DGAlgCohRing(C::SimplicialCochainComplex)

Given a cochain complex `C` with a structure as a differential graded algebra 
(via implementation of the internal function `mul_cochains`) this creates the 
associated cohomology ring.
"""
  function DGAlgCohRing(C::SimplicialCochainComplex)
    T = elem_type(base_ring(C))
    return new{T}(C)
  end
end

@doc raw"""
    mul_cochains(C::AbsHyperComplex, a, p::Int, b, q::Int)

Cochain complexes to be used in `DGAlgCohRing`s need to be endowed with 
the structure of a DG-algebra. To achieve this, the programmer needs to overwrite 
the method for this function in accordance with `https://stacks.math.columbia.edu/tag/061V`.

In short this must implement the multiplication of two representatives `a` and `b` of 
cohomology classes in `C` of degrees `p` and `q` and return a representative of a 
class in degree `p + q`. 
"""
function mul_cochains(C::AbsHyperComplex, a, p::Int, b, q::Int)
  error("not implemented for complexes of type $(typeof(C)); see the source code for more information")
end

@doc raw"""
    simplicial_cochain_complex(A::DGAlgCohRing)

Return the internal cocomplex (which needs to implement a structure as a DG-algebra).
"""
simplicial_cochain_complex(A::DGAlgCohRing) = A.C


@doc raw"""
    graded_parts(A::DGAlgCohRing)

Return a list of the cohomology groups for `A`.
"""
function graded_parts(A::DGAlgCohRing)
  if !isdefined(A, :graded_parts)
    list = [homology(A.C, 0)[1]]
    i = 1
    while can_compute_index(A.C, (i,))
      push!(list, homology(A.C, i)[1])
      i += 1
    end
    A.graded_parts = list
  end
  return A.graded_parts
end

@doc raw"""
    graded_part(A::DGAlgCohRing, i::Int)

Return the `i`-th cohomology group of the cocomplex for `A`.
"""
function graded_part(A::DGAlgCohRing, i::Int)
  return graded_parts(A)[i+1]
end

mutable struct DGAlgCohRingElem{T} <: NCRingElem
  parent::DGAlgCohRing{T}
  # Elements can be represented in three ways:
  #  1) no fields assigned (the zero element)
  #  2) `homog_elem` and `homog_deg` assigned (homogeneously stored element)
  #  3) a dictionary `p => v` where the degree points to the cohomology 
  #     class for the degree `p` part of this element.
  # The values of `coeff` must not be zero. The zero element in 
  # the ring must be stored as as the empty element. 
  #
  # The distinction whether or not an element is stored as a homogeneous 
  # element can be made by whether or not the field `.homog_elem` is set 
  # to `nothing`. The field `.coeff` is not assigned if it's not needed. 
  homog_elem::Union{SubquoModuleElem{T}, Nothing}
  homog_deg::Union{Int, Nothing}
  coeff::Dict{Int, SubquoModuleElem{T}}

  # Constructor for homogeneous elements
  function DGAlgCohRingElem(
      A::DGAlgCohRing{T}, 
      p::Int, v::SubquoModuleElem{T}
    ) where {T}
    @assert parent(v) === graded_part(A, p)
    return new{T}(A, v, p)
  end

  # Constructor for mixed elements
  function DGAlgCohRingElem(
      A::DGAlgCohRing{T}, 
      coeff::Dict{Int, SubquoModuleElem{T}};
      check::Bool=true
    ) where {T}
    @check all(parent(v) === graded_part(A, p) for (p, v) in coeff) "degree/element incompatibility"
    return new{T}(A, nothing, nothing, coeff)
  end

  # Constructor for the zero
  function DGAlgCohRingElem(
      A::DGAlgCohRing{T}
    ) where {T}
    return new{T}(A, nothing, nothing)
  end
end

parent(a::DGAlgCohRingElem) = a.parent
zero(A::DGAlgCohRing) = DGAlgCohRingElem(A)
zero(a::DGAlgCohRingElem) = zero(parent(a))
(A::DGAlgCohRing)() = zero(A)

function one(A::DGAlgCohRing)
  DGAlgCohRingElem(A, 0, sum(gens(homology(simplicial_cochain_complex(A), 0)[1])))
end


@doc raw"""
    graded_part(a::DGAlgCohRingElem, p::Int)

Return the homogeneous component of `a` of cohomological degree `p`.
"""
function graded_part(a::DGAlgCohRingElem, p::Int)
  A = parent(a)
  is_zero(a) && return zero(graded_part(A, p))
  if !isnothing(a.homog_elem)
    p == a.homog_deg || return zero(graded_part(A, p))
    return a.homog_elem
  end
  return get(a.coeff, p, zero(graded_part(A, p)))
end

# Coercion from base ring
function (A::DGAlgCohRing{T})(c::T) where T
  one_A = one(A)
  one_A.homog_elem = c*one_A.homog_elem
  return one_A
end

# General coercion
function (A::DGAlgCohRing)(c)
  R = base_ring(A)
  return(A(R(c)))
end

function (A::DGAlgCohRing)(c::DGAlgCohRingElem)
  parent(c) === A || error("wrong parent")
  return(c)    
end

function is_zero(a::DGAlgCohRingElem)
  if isnothing(a.homog_elem)
    !isdefined(a, :coeff) && return true
    isempty(a.coeff) && return true
    return all(iszero(b) for (_, b) in a.coeff)
  end
  return is_zero(a.homog_elem)
end

function deepcopy_internal(a::DGAlgCohRingElem, d::IdDict)
  result = parent(a)()
  if !isnothing(a.homog_elem)
    result.homog_elem = deepcopy_internal(a.homog_elem, d)
    result.homog_deg = deepcopy(a.homog_deg)
    return result
  end
  if isdefined(a, :coeff)
    result.coeff = deepcopy_internal(a.coeff, d)
  end
  return result
end

function +(a::DGAlgCohRingElem{T}, b::DGAlgCohRingElem{T}) where {T}
  return add!(deepcopy(a), b)
end

function add!(a::DGAlgCohRingElem{T}, b::DGAlgCohRingElem{T}) where {T}
  @assert parent(a) === parent(b) "parent mismatch"
  A = parent(a)
  is_zero(a) && return deepcopy(b)
  is_zero(b) && return a
  if isnothing(a.homog_elem)
    !isdefined(a, :coeff) && return deepcopy(b) # a was the zero element in the end
    if isnothing(b.homog_elem)
      for (q, v) in b.coeff
        w = get(a.coeff, q, nothing)
        if isnothing(w)
          a.coeff[q] = v
        else
          res = w + v
          if is_zero(res)
            delete!(a.coeff, q)
          else
            a.coeff[q] = res
          end
        end
      end
      isempty(a.coeff) && return zero(A)
      return a
    else # b is homogeneous
      q = b.homog_deg
      w = get(a.coeff, q, nothing)
      if isnothing(w)
        a.coeff[q] = b.homog_elem
      else
        res = w + b.homog_elem
        if is_zero(res)
          delete!(a.coeff, q)
        else
          a.coeff[q] = res
        end
      end
      isempty(a.coeff) && return zero(A)
      return a
    end
  else # a is homogeneous
    q = a.homog_deg
    result = zero(A) # we allocate new because we can not delete a.homog_elem
    if isnothing(b.homog_elem) # b is not homogeneous
      result.coeff = deepcopy(b.coeff) # result will not be homogeneously stored
      w = get(result.coeff, q, nothing)
      if isnothing(w)
        result.coeff[q] = a.homog_elem
      else
        res = w + a.homog_elem
        if iszero(res)
          delete!(result.coeff, q)
        else
          result.coeff[q] = res
        end
      end
      isempty(result.coeff) && return zero(A)
      return result
    else # b is homogeneous
      p = b.homog_deg
      if p == q
        result.homog_deg = p
        result.homog_elem = a.homog_elem + b.homog_elem
        is_zero(result.homog_elem) && return zero(A)
        return result
      else
        result.coeff = Dict{Int, SubquoModuleElem{T}}([q => a.homog_elem, p => b.homog_elem])
        isempty(result.coeff) && return zero(A)
        return result
      end
    end
  end
end

# get i-th (module) generator of H^p
Base.getindex(A::DGAlgCohRing, p::Int, i::Int) = DGAlgCohRingElem(A, p, gens(homology(A.C, p)[1])[i])

is_homogeneous_normalized(a::DGAlgCohRingElem) = !isnothing(a.homog_elem)

is_homogeneous_denormalized(a::DGAlgCohRingElem) = isnothing(a.homog_elem) && isone(length(a.coeff))

is_homogeneous(a::DGAlgCohRingElem) = iszero(a) || is_homogeneous_normalized(a) || is_homogeneous_denormalized(a)

function degree(a::DGAlgCohRingElem)
  if is_homogeneous_normalized(a)
    return a.homog_deg
  elseif is_homogeneous_denormalized(a)
    return only(keys(a.coeff))
  else
    throw(ArgumentError("element is not homogeneous"))
  end
end

-(a::DGAlgCohRingElem) = is_homogeneous_normalized(a) ? DGAlgCohRingElem(parent(a), a.homog_deg, -a.homog_elem) : sum(-b for b in homogeneous_parts(a); init=zero(parent(a)))

-(a::DGAlgCohRingElem, b::DGAlgCohRingElem) = a + (-b)

# extract homogeneous parts. We do a set, because we need it for ==
function homogeneous_parts(a::DGAlgCohRingElem)
  A = parent(a)
  is_zero(a) && return Set{elem_type(graded_part(A, 0))}()
  isnothing(a.homog_elem) && return Set(DGAlgCohRingElem(parent(a), i, m) for (i,m) in a.coeff)
  return Set([a])
end

# distribute over homogeneous parts
function *(a::DGAlgCohRingElem, b::DGAlgCohRingElem)
  @req parent(a) === parent(b) "parent mismatch"
  if is_zero(a) || is_zero(b)
    return zero(a)
  end
  if !is_homogeneous_normalized(a)
    return sum(ah*b for ah in homogeneous_parts(a); init=zero(a))
  end
  # a is homogeneous
  if !is_homogeneous_normalized(b)
    return sum(a*bh for bh in homogeneous_parts(b); init=zero(a))
  end

  # both a and b are homogeneous
  return mul_homog(a, b)
end

# Multiplication on homogeneous parts
function mul_homog(a, b)
  @req parent(a) === parent(b) "parent mismatch"
  p = a.homog_deg
  q = b.homog_deg
  A = parent(a)
  C = simplicial_cochain_complex(A)
  !can_compute_index(C, p+q) && return zero(A)
  H = homology(C, p+q)[1]
  cochain = mul_cochains(C, repres(a.homog_elem), p, repres(b.homog_elem), q)
  return DGAlgCohRingElem(A, p+q, H(cochain; check=false))
end

# Parent and element types
elem_type(::Type{DGAlgCohRing{T}}) where {T} = DGAlgCohRingElem{T}
parent_type(::Type{DGAlgCohRingElem{T}}) where {T} = DGAlgCohRing{T}


# Base ring and base ring type
base_ring(C::DGAlgCohRing) = base_ring(simplicial_cochain_complex(C))
base_ring_type(::Type{DGAlgCohRingElem{T}}) where {T} = parent_type(T)
base_ring_type(::Type{DGAlgCohRing{T}}) where {T} = parent_type(T)
# Equality for homogeneous elements is straight foward; for inhomogeneous, do it by sets of homogeneous parts
function ==(a::DGAlgCohRingElem, b::DGAlgCohRingElem)
  if is_homogeneous_normalized(a) && is_homogeneous_normalized(b)
    return a.homog_elem == b.homog_elem && a.homog_deg == b.homog_deg
  else
    return homogeneous_parts(a) == homogeneous_parts(b)
  end
end

function Base.hash(a::DGAlgCohRingElem, u::UInt)
  A = parent(a)
  r = 0x60125ab8e8cd44ca
  if isdefined(a, :homog_elem)
    r = xor(r, hash(a.homog_deg, u))
    r = xor(r, hash(a.homog_elem, u))
    return r
  end

  !isdefined(a, :coeff) && return xor(r, u)
  
  # we need to catch the case where homogeneous elements use the dictionary
  if isone(length(a.coeff))
    p, c = only(a.coeff)
    r = xor(r, hash(p, u))
    r = xor(r, hash(c))
    return r
  end

  return hash(a.coeff, u)
end

# This only says that not every cohomology ring is a domain, it does not look at a specific ring.
function is_domain_type(::Type{DGAlgCohRingElem})
  return false
end

function is_exact_type(::Type{DGAlgCohRingElem{T}}) where T
  return is_exact_type(T)
end

Base.show(io::IO, a::DGAlgCohRingElem) = Base.show(io, MIME"text/plain"(), a)

# Show cohomology ring element as string.
# If context provides :parens=>true, wrap result in parenthesis
function Base.show(io::IO, mime::MIME"text/plain", a::DGAlgCohRingElem)
  if get(io, :parens, false)
    print(io, "(")
  end
  if iszero(a)
    print(io, "0")
  elseif is_homogeneous_normalized(a)
    # defer printing on the actual cochain to the cochain complex
    show_elem(io, simplicial_cochain_complex(parent(a)), repres(a.homog_elem), a.homog_deg)
    print(io, " + ", a.homog_deg, "-coboundaries")
  elseif is_homogeneous_denormalized(a)
    print(IOContext(io, :parens=>false), only(homogeneous_parts(a)))
  else
    join(IOContext(io, :parens=>true), homogeneous_parts(a), " + ")
  end
  if get(io, :parens, false)
    print(io, ")")
  end
end

# The canonical unit is the one of the cohomology ring. 
canonical_unit(A::DGAlgCohRing) = one(A)

# isone checks if the element is the natural one element
isone(a::DGAlgCohRingElem) = (a == one(parent(a)))

# Interface specific to noncommutative rings
function divexact_left(a::DGAlgCohRingElem, b::DGAlgCohRingElem)
  if a == b 
    return one(parent(a))
  end 
  error("no exact quotient exists")
end

function divexact_right(b::DGAlgCohRingElem, a::DGAlgCohRingElem)
  if a == b 
    return one(parent(a))
  end 
  error("no exact quotient exists")
end


# we need exponentiation for tests
function ^(a::DGAlgCohRingElem, n::Int)
  n >= 0 || error("negative exponent")
  is_one(n) && return deepcopy(a)
  is_zero(n) && return one(a)
  p, r = divrem(n, 2)
  return a^p * a^(p+r)
end

# needed for the test suite but not to be used elsewhere!
function generate_homogeneous_element(R::Oscar.DGAlgCohRing{ZZRingElem})
  degree = rand(1:dim(simplicial_complex(simplicial_cochain_complex(R)))+1) # indexing starts at 1 (for degree 0)
  n_gens = rand(0:2*length(gens(Oscar.graded_parts(R)[degree])))
  x = zero(R)
  for i=1:n_gens
    n = rand(-100:100)
    x = x+n*R[degree-1,rand(1:length(gens(Oscar.graded_parts(R)[degree])))]
  end
  return x
end

@doc raw"""
    set_volume_form!(A::DGAlgCohRing, v::DGAlgCohRingElem)
    set_volume_form!(A::DGAlgCohRing)

A volume form `v` is a cohomology class in the top-dimensional non-trivial degree of a 
cohomology ring which generates that degree as a free module of rank one over the 
`coefficient_ring`. 

In case no element `v` is provided, a heuristic is applied to determine a random 
suitable element. This heuristic may fail depending on existence of implementations 
for `ModuleFP`s for the chosen `coefficient_ring`. 
"""
function set_volume_form!(A::DGAlgCohRing, v::DGAlgCohRingElem)
  isdefined(A, :volume_form) && error("volume form is already defined; changes are not possible")
  @assert is_homogeneous(v)
  d = degree(v)
  A.volume_form = (d, graded_part(v, d))
  return A
end

function set_volume_form!(A::DGAlgCohRing)
  @assert has_upper_bound(A.C) "upper bound for the underlying complex needs to be known"
  b = upper_bound(A.C)
  while (is_zero(A.C[b]) || is_zero(homology(A.C, b)[1])) && b >= 0
    b -= 1
  end
  b == 0 && error("no non-zero cohomology group found")
  g = small_generating_set(A, b)
  @assert is_one(length(g)) "top cohomology group could not be identified as free of rank one"
  set_volume_form!(A, only(g))
end

@doc raw"""
    volume_form(A::DGAlgCohRing)

If the top-dimensional, non-trivial graded part of `A` is freely generated 
of rank one, one may choose a `volume_form`, i.e. a generator of that graded 
part over the `coefficient_ring` via `set_volume_form!`. This returns the once 
chosen form if it was set.
"""
function volume_form(A::DGAlgCohRing)
  !isdefined(A, :volume_form) && error("volume form has not been chosen; use `set_volume_form!` for that")
  d, v = A.volume_form
  DGAlgCohRingElem(A, d, v)
end

@doc raw"""
    integral(a::DGAlgCohRingElem)

Assuming the `parent` of `a` is provided with a chosen `volume_form`, return 
the integral of `a` with respect to that form.

Use `set_volume_form!` to choose a volume form on a `DGAlgCohRing`. 
"""
function integral(a::DGAlgCohRingElem)
  A = parent(a)
  inc = volume_form_inclusion(A)
  d, _ = A.volume_form # grabbing the field is OK as the previous call makes sure it's set
  return preimage(inc, graded_part(a, d))[1]
end

function volume_form_inclusion(A::DGAlgCohRing)
  if !isdefined(A, :vol_form_inc)
    vol = volume_form(A)
    d = degree(vol)
    v = graded_part(vol, d)
    _, inc = sub(graded_part(A, d), [v])
    A.vol_form_inc = inc
  end
  return A.vol_form_inc
end

function small_generating_set(A::DGAlgCohRing{T}, d::Int) where T
  if !isdefined(A, :small_gens)
    A.small_gens = Dict{Int, Vector{SubquoModuleElem{T}}}()
  end
  g = get!(A.small_gens, d) do
    Hd = graded_part(A, d)
    M, iso = simplify(Hd)
    return iso.(gens(M))
  end
  return elem_type(A)[DGAlgCohRingElem(A, d, v) for v in g]
end


