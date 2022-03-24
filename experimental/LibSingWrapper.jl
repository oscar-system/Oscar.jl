import Singular.LibSing: *

export milnor, tjurina, global_milnor_number

### The Milnor algebra of an isolated hypersurface singularity f

@Markdown.doc """
    milnor(f::MPolyElem; point::Vector=[])

For an isolated hypersurface singularity `f` at a point `point`, 
compute the triple `(gb, μ, g)` where `gb` is a standard basis for 
the jacobian ideal, `μ` is the Milnor number, and `g` is a ``𝕜``-base 
for the Milnor algebra.

The default argument for `point` is the origin.
"""
function milnor(f::MPolyElem; point::Vector=[])
  iszero(f) && error("polynomial is zero")
  R = parent(f)
  kk = coefficient_ring(f)
  # check whether we are working over a field of coefficients
  typeof(kk)<:AbstractAlgebra.Field || error("coefficient domain is not a field")

  # We can not provide a full default vector of zeroes in the 
  # signature of the function, because we would need a full 
  # instance of the field at that point. Hence, we have to 
  # do that at a later point here.
  k_point = kk.(point)
  if length(point) != 0
    # check for the right number of coordinates
    length(k_point) == ngens(R) || error("vector does not have the correct length")
  else
    # fill up the default otherwise
    k_point = elem_type(kk)[zero(kk) for i in 1:ngens(R)]
  end

  # Set up the localization at `point` so that we have 
  # the local orderings at our disposal
  W = Localization(MPolyComplementOfKPointIdeal(R, k_point))

  # compute the jacobian ideal
  Df = jacobi_matrix(f)::AbstractAlgebra.Generic.MatSpaceElem{typeof(f)}
  I = ideal(R, [Df[i, 1] for i in 1:nrows(Df)])
  
  # transfer to the localized ring since otherwise we could not work with 
  # local orderings
  WI = W(I)
  # compute a standard basis w.r.t a local ordering 
  stdI = groebner_basis(WI)

  # carry out the computations on the Singular side
  Singular.dimension(singular_gens(stdI)) == 0 || error("singularity is not isolated")
  milnor_number = Singular.vdim(singular_gens(stdI)) 
  # make the conversion to the Oscar side 
  milnor_kbase = [to_oscar_side(stdI, g) for g in gens(Singular.kbase(singular_gens(stdI)))]
  return (stdI, milnor_number, milnor_kbase)
end

### The Tjurina number  of an isolated hypersurface singularity f
@Markdown.doc """
    tjurina(f::MPolyElem; point::Vector=[])

For an isolated hypersurface singularity `f` at a point `point`,
compute the triple `(gb,\tau,g)` where `gb` is a standard basis for the
Tjurina ideal. `\tau` is the Tjurina number, and `g` is a ``k``-base of the
Tjurina algebra.

The default argument for `point` is the origin.
"""
function tjurina(f::MPolyElem; point::Vector=[])
  iszero(f) && error("input polynomial is zero")
  R = parent(f)
  kk = coefficient_ring(f)
  # make sure that kk is a field
  typeof(kk)<:AbstractAlgebra.Field || error("coefficient domain is not a field")

  # Set sanity check on point or creation of vector of point, if default point 
  k_point = kk.(point)
  if length(point) != 0
     # right number of coordinates?
     length(k_point) == ngens(R) || error("vector does not have correct length")
  else
     # default case, need to set up point for origin
     k_point = elem_type(kk)[zero(kk) for i in 1:ngens(R)]
  end

  # localize at point (setting up stage for Singular)
  RL = Localization(MPolyComplementOfKPointIdeal(R,k_point))  

  # compute the jacobian and tjurina ideals
  Jf = jacobi_matrix(f)::AbstractAlgebra.Generic.MatSpaceElem{typeof(f)}
  I=ideal(R,[Jf[i,1] for i in 1:nrows(Jf)]) + ideal(R,[f])
  # now move to localization RL of R and do standard basis there
  IL=RL(I)
  stdIL = groebner_basis(IL)
  
  # now fill in the return values
  Singular.dimension(singular_gens(stdIL)) == 0 || error("singularity is not isolated")
  tjurina_number = Singular.vdim(singular_gens(stdIL)) 
  # make the conversion to the Oscar side 
  tjurina_kbase = [to_oscar_side(stdIL, g) for g in gens(Singular.kbase(singular_gens(stdIL)))]
  return (stdIL, tjurina_number, tjurina_kbase)
end

### The global Milnor number of a hypersurface with at most isolated singularities
@Markdown.doc """
    global_milnor_number(f::MPolyElem)

Computes the global Milnor number of the affine hypersurface 
defined by `f` with at most isolated singularities.
"""
function global_milnor_number(f::MPolyElem)
  iszero(f) && error("polynomial is zero")
  R = parent(f)
  kk = coefficient_ring(f)
  # check whether we are working over a field of coefficients
  typeof(kk)<:AbstractAlgebra.Field || error("coefficient domain is not a field")

  # compute groebner basis of the jacobian ideal
  Df = jacobi_matrix(f)::AbstractAlgebra.Generic.MatSpaceElem{typeof(f)}
  #I = ideal(R, [Df[i, 1] for i in 1:nrows(Df)])

  RS=Oscar.singular_ring(R)
  # I needs to be passed to RS where we can make a Singular Standard basis..
  Ising = Singular.Ideal(RS, [RS(Df[i, 1]) for i in 1:nrows(Df)])
  stdI = Singular.std(Ising)

  # set up the resulting data
  Singular.dimension(stdI) == 0 || error("non-isolated singularities detected")
  milnor_number = Singular.vdim(stdI) 
  return (milnor_number)
end

@Markdown.doc """
    milnor(f::Vector{T}) where {T<:MPolyElem}

For an isolated complete intersection singularity given by 
the regular sequence ``f₁,…,fₖ`` compute the Milnor number 
at `point` by means of the Le-Greuel-formula.
"""
function milnor(f::Vector{T}; point::Vector=[]) where {T<:MPolyElem}
  length(f) == 0 && error("not an ICIS")
  R = parent(f[1])
  k = length(f)
  for i in 2:k
    parent(f[i]) == R || error("elements do not have the same parent")
  end
  n = ngens(R)

  kk = coefficient_ring(R)
  # check whether we are working of a field of coefficients
  typeof(kk)<:AbstractAlgebra.Field || error("coefficient domain is not a field")

  # We can not provide a full default vector of zeroes in the 
  # signature of the function, because we would need a full 
  # instance of the field at that point. Hence, we have to 
  # do that at a later point here.
  k_point = kk.(point)
  if length(point) != 0
    # check for the right number of coordinates
    length(k_point) == ngens(R) || error("vector does not have the correct length")
  else
    # fill up the default otherwise
    k_point = elem_type(kk)[zero(kk) for i in 1:ngens(R)]
  end

  rem_f = f
  summands = Int[]
  W = Localization(MPolyComplementOfKPointIdeal(R, k_point))
  sign = 1
  for i in k:-1:1
    Df_rem = jacobi_matrix(rem_f)
    Df_min = minors(Df_rem, i)
    pop!(rem_f)
    Df_min_ext = vcat(Df_min, rem_f)
    I = ideal(R, Df_min_ext)
    WI = W(I)
    stdI = groebner_basis(WI)
    Singular.dimension(singular_gens(stdI)) == 0 || error("not an ICIS or not in general position")
    b = Singular.vdim(singular_gens(stdI))
    #b >= 0 || error("not an ICIS or not in general position")
    push!(summands, sign*b)
    sign = -sign
  end
  return sum(summands)
end
