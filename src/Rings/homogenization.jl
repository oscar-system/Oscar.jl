export dehomogenization, homogenization

###############################################################################
# Homogenization and dehomogenization
################################################################################
# dehomogenization with respect to the specified variable into the specified
# ring

@doc Markdown.doc"""
    dehomogenization([r::MPolyRing], F::Oscar.MPolyElem_dec, i::Int)

Return the polynomial obtained by dehomogenization of `F` with respect to the `i`th variable of `parent(F)`. The new polynomial is in `r` if specified, or in a new ring with one variable less otherwise.
"""
function dehomogenization(r::MPolyRing{S}, F::Oscar.MPolyElem_dec{S}, i::Int) where S <: FieldElem
  @assert ishomogenous(F)
  R = parent(F)
  nvars(R) -1 == nvars(r) || error("Incompatible number of variables")
  V = gens(r)
  insert!(V, i, r(1))
  phi = hom(R.R, r, V)
  return phi(F.f)
end

################################################################################
# dehomogenization without specifying the ring with repect to the specified variable

function dehomogenization(F::Oscar.MPolyElem_dec, i::Int)
  R = parent(F)
  A = String.(symbols(R))
  r = PolynomialRing(R.R.base_ring, deleteat!(A, i))
  dehomogenization(r[1], F, i)
end

################################################################################
# non decorated version

function dehomogenization(F::Oscar.MPolyElem, i::Int)
  R = parent(F)
  A = grade(R)
  dehomogenization(A(F), i)
end

################################################################################

function dehomogenization(r::MPolyRing{S}, F::Oscar.MPolyElem{S}, i::Int) where S <: FieldElem
  R = parent(F)
  A = grade(R)
  dehomogenization(r, A(F), i)
end

################################################################################
# homogenization

 @doc Markdown.doc"""
     homogenization(R::MPolyRing{S}, F::MPolyElem{S}, i::Int) where S <: FieldElem

Return the homogenization of the polynomial `F` in `R` using the `i`th variable
of `R`, and mapping the variables of `parent(F)` to the variables of `R`
respecting the same order.
 """
 function homogenization(R::MPolyRing{S}, F::Oscar.MPolyElem{S}, i::Int) where S <: FieldElem
   r = parent(F)
   V = gens(R)
   W = [V[j]//V[i] for j=1:nvars(R)]
   deleteat!(V, i)
   phi = hom(r, R, V)
   G = phi(F)
   d = total_degree(F)
   P = R[i]^d*evaluate(G, W)
   return numerator(P)
 end

################################################################################

@doc Markdown.doc"""
      homogenization(F::Oscar.MPolyElem; variable::String="x0", position::Bool = true)

Return the homogenization of the polynomial `F` in a ring with one additional
variable (named `x0` if not specified). If no position is specified, the additional
variable is the first variable of the new ring, if `position = false` is specified,
the new variable is the last variable of the new ring.
 """
 function homogenization(F::Oscar.MPolyElem; variable::String="x0", position::Bool = true)
   r = parent(F)
   A = String.(symbols(r))
   if position == false
     push!(A, variable)
     i = length(A)
   else
     prepend!(A, [variable])
     i = 1
   end
   R = PolynomialRing(r.base_ring, A)
   homogenization(R[1], F, i)
 end

################################################################################

@doc Markdown.doc"""
    homogenization(R::MPolyRing{S}, V::Vector{T}, i::Int) where {S, T <: MPolyElem{S}}

Return the homogenization of the elements of `V` in `R` using the `i`th variable of `R`.
"""
function homogenization(R::MPolyRing{S}, V::Vector{T}, i::Int) where {S, T <: MPolyElem{S}}
  return [homogenization(R, v, i) for v in V]
end

################################################################################

@doc Markdown.doc"""
      homogenization(V::Vector{S}; variable::String = "x0", position::Bool = true) where S <: MPolyElem

Return the homogenization of the elements of `V` in a ring with one additional
variable (named `x0` if not specified). If no position is specified, the additional
variable is the first variable of the new ring, if `position = false` is specified,
the new variable is the last variable of the new ring.
 """
function homogenization(V::Vector{S}; variable::String = "x0", position::Bool = true) where S <: MPolyElem
  @assert !isempty(V)
  r = parent(V[1])
  A = String.(symbols(r))
  if position == false
    push!(A, variable)
    i = length(A)
  else
    prepend!(A, [variable])
    i = 1
  end
  R = PolynomialRing(r.base_ring, A)
  return homogenization(R[1], V, i)
end

################################################################################
@doc Markdown.doc"""
    homogenization(R::MPolyRing{S}, I::MPolyIdeal{T}, i::Int, ordering::Symbol = :degrevlex) where {S, T <: MPolyElem{S}}

Return the homogenization of the ideal `I` in the ring `R`, using the ordering
degrevlex if not specified otherwise and using the `i`th variable of `R`.
"""
function homogenization(R::MPolyRing{S}, I::MPolyIdeal{T}, i::Int, ordering::Symbol = :degrevlex) where {S, T <: MPolyElem{S}}
  return ideal(R, homogenization(R, groebner_basis(I, ord = ordering), i))
end

################################################################################

@doc Markdown.doc"""
      homogenization(F::Oscar.MPolyElem; variable::String="x0", position::Bool = true, ordering::Symbol = :degrevlex)

Return the homogenization of the ideal `I`, using the ordering degrevlex if not
specified otherwise, in a ring with one additional variable (named `x0` if not
specified). If no position is specified, the additional variable is the
first variable of the new ring, if `position = false` is specified, the new
variable is the last variable of the new ring.
 """
function homogenization(I::MPolyIdeal; variable::String = "x0", position::Bool = true, ordering::Symbol = :degrevlex)
  r = base_ring(I)
  A = String.(symbols(r))
  if position == false
    push!(A, variable)
    i = length(A)
  else
    prepend!(A, [variable])
    i = 1
  end
  R = PolynomialRing(r.base_ring, A)
  return homogenization(R[1], I, i, ordering)
end
