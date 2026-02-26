# Ad hoc documentation (just as a comment currently)

# (BUGGY) PROTOTYPE for generating random values

# Summary of features for this prototype (in no special order)
# (1) explicit distribution objects
# (2) RAND(parent, distrib)  or  RAND(RNG, parent, distrib)
# (3) DistrList(parent) gives list of possible distribution objects for parent
# (4) RAND_(parent)  uses  DefaultDistribution(parent)
# (5) WidenDistr(D)  produces a possibly "wider" distribution
# (6) ApproxCard  gives approximate cardinality of range of possible random values

# Points (4), (5) and (6) are still "dodgy".
# Function/struct and other names are _ad hoc_; similarly for fn prototypes.

#= EXAMPLES

RAND(ZZ, UniformIntegerDistr(-9, 9))

QQx,x = QQ["x"]
DistrList(QQx)
SmallCoeffDistr = RationalDistr(UniformIntegerDistr(-8,8), UniformIntegerDistr(1,4))
SmallPolyDistr = UnivariatePolyDistr(3, SmallCoeffDistr)
LargeCoeffDistr = RationalDistr(UniformIntegerDistr(-999,999), UniformIntegerDistr(1,999))
LargePolyDistr = UnivariatePolyDistr(8, LargeCoeffDistr)
[ RAND(QQx, SmallPolyDistr)  for i in 1:5 ]
[ RAND(QQx, LargePolyDistr)  for i in 1:5 ]

# ---------------------------------
# "lazy" function RAND_
# YOU PROBABLY DO NOT WANT TO USE RAND_

[ RAND_(ZZ)  for i in 1:6 ]

QQxy,(x,y) = QQ["x", "y"]
[ RAND_(QQxy)  for i in 1:6 ]
=#

# ==================================================================
# Some basic "distribution" object types

abstract type AbstractRndDistr end

function ContainsZero(D::AbstractRndDistr) end

function ApproxCard(D::AbstractRndDistr)
  return D.ApproxCardinality
end

function PostFilter(predicate, RNG::Random.AbstractRNG, parent, distr::AbstractRndDistr)
  # Generate a random value satisfying predicate: give up after 50 attempts if no suitable value found.
  # Assuming we have at least 50% chance of generating a suitable value, this will fail needlessly
  # only very rarely.
  count_trials = 0
  while true
    val = RAND(RNG, parent, distr)
    if predicate(val)
      return val
    end
    count_trials += 1
    (count_trials > 50) && error("Failed to find random element satsfying given predicate")
  end
  # never get here
end

# THIS DID NOT WORK AS HOPED: e.g. try on identity_matrix(QQ,2)
# # DODGY GENERAL FUNCTION
# function DistrList(::T) where {T <: Set}
#   return DistrList(T)
# end

# ------------------------------------------------------------------
# ----- UniformDistr -----

struct UniformDistr <: AbstractRndDistr
  ApproxCardinality::Int  # UH OH  !!! BIG PROBLEM !!!

  function UniformDistr()
    return new()
  end
end

function ContainsZero(D::UniformDistr)
  return true # assuming it is applied to a finite ring/field
end

function RAND(RNG::Random.AbstractRNG, parent::T, D::UniformDistr) where { T <: FinField }
  return parent(rand(RNG, 0:characteristic(parent)-1))  # !!! WRONG !!! (valid only for prime subfield)
end
  
function RAND(RNG::Random.AbstractRNG, parent::T, D::UniformDistr) where { T <: zzModRing }
  return parent(rand(RNG, 0:characteristic(parent)-1))  # !!! WRONG !!! (valid on for base subring)
end
  
function RAND(RNG::Random.AbstractRNG, parent::T, D::UniformDistr) where { T <: ZZModRing }
  return parent(rand(RNG, 0:characteristic(parent)-1))  # !!! WRONG !!! (valid on for base subring)
end
  
function DistrList(::T) where { T <: FinField }
  return [UniformDistr];
end

function DistrList(::T) where { T <: zzModRing }
  return [UniformDistr];
end

function DistrList(::T) where { T <: ZZModRing }
  return [UniformDistr];
end


# ------------------------------------------------------------------
# ----- UniformIntegerDistr -----

struct UniformIntegerDistr <: AbstractRndDistr
  lwb::ZZRingElem
  upb::ZZRingElem # always have lwb <= upb
  ApproxCardinality::Int

  function UniformIntegerDistr(lo::ZZRingElem, hi::ZZRingElem)
    (lo <= hi) || error("Empty range")
    return new(lo,hi, hi-lo+1)
  end
end

# Constructor from machine integers:
function UniformIntegerDistr(lo::Int, hi::Int)
  return UniformIntegerDistr(ZZ(lo), ZZ(hi))
end

function ContainsZero(D::UniformIntegerDistr)
  return (D.lwb <= 0 && D.upb >= 0)
end

function width(D::UniformIntegerDistr)
  return D.upb-D.lwb
end

function RAND(RNG::Random.AbstractRNG, parent, D::UniformIntegerDistr)
  # Could behave "oddly" e.g. with input parent=ZZ/(101) and UnifRange=[0,101], so zero has a double probability!
  return parent(rand(RNG, (D.lwb):(D.upb)))
end
  
function DistrList(::ZZRing)
  return [UniformIntegerDistr];
end

# ------------------------------------------------------------------
# ----- RationalDistr -----

struct RationalDistr <: AbstractRndDistr
  NumerDistr::UniformIntegerDistr
  DenomDistr::UniformIntegerDistr  # range is always positive
  ApproxCardinality::Int64

  function RationalDistr(numer::UniformIntegerDistr, denom::UniformIntegerDistr)
    (denom.lwb > 0) || error("Denom range must be positive")
    # !!! IMPROVE LINE BELOW !!!
    card = div(ZZ(3)*ApproxCard(numer)*ApproxCard(denom), 5)
    card = Int64(clamp(card,0,typemax(Int64)))
    return new(numer, denom, card)
  end
end

function ContainsZero(D::RationalDistr)
  return ContainsZero(D.NumerDistr)
end

function RAND(RNG::Random.AbstractRNG, parent, distr::RationalDistr)
  N = RAND(RNG, parent, distr.NumerDistr)
  D = RAND(RNG, parent, distr.DenomDistr)
  return N/D # could give error (div by zero)
end

function DistrList(::QQField)
  return [RationalDistr];
end

# ------------------------------------------------------------------
# ----- (dense) UnivariatePoly -----

struct UnivariatePolyDistr <: AbstractRndDistr
  deg::Int  # *exact* degree of result, d >= 0
  CoeffDistr::AbstractRndDistr
  ApproxCardinality::Int64

  function UnivariatePolyDistr(d::Int, coeff::AbstractRndDistr)
    (d >= 0) || error("Degree must be non-negative")
    CardLC = ApproxCard(coeff); if ContainsZero(coeff)  CardLC -= 1; end
    (CardLC > 0) || error("Coeff distribution must allow non-zero values")
    card = clamp(CardLC*ZZ(ApproxCard(coeff))^d, 0, typemax(Int64))
    card = Int64(card)
    return new(d, coeff, card)
  end
end

function ContainsZero(D::UnivariatePolyDistr)
  return false
end

function RAND(RNG::Random.AbstractRNG, parent::T, distr::UnivariatePolyDistr) where { T <: PolyRing }
  d = distr.deg
  x = gen(parent)
  k = coefficient_ring(parent)
  LC = PostFilter(!is_zero, RNG, k, distr.CoeffDistr)
  poly = zero(parent)
  for i in 0:d-1
    poly += (x^i) * RAND(RNG, k, distr.CoeffDistr)
  end
  return LC*x^d + poly
end

function DistrList(::T) where { T <: PolyRing }
  return [UnivariatePolyDistr];
end

# ------------------------------------------------------------------
# ----- MultivariatePoly -----

struct MultivariatePolyCubeDistr <: AbstractRndDistr
  deg::Int  # each variable appears to degree at most deg
  nterms::Int # polynomial will have at most nterms terms
  CoeffDistr::AbstractRndDistr
  ApproxCardinality::Int64

  function MultivariatePolyCubeDistr(d::Int, t::Int, coeff::AbstractRndDistr)
    (d >= 0) || error("Degree must be non-negative")
    (t >= 0) || error("Number of terms must be non-negative")
    CardLC = ApproxCard(coeff); if ContainsZero(coeff)  CardLC -= 1; end
    (CardLC > 0) || error("Coeff distribution must allow non-zero values")
    card = clamp(ZZ(CardLC)^t, 0, typemax(Int64))  # wasteful power 
    card = Int64(card)
    return new(d, t, coeff, card)
  end
end

function ContainsZero(D::MultivariatePolyCubeDistr)
  return ContainsZero(D.CoeffDistr)
end

function RAND(RNG::Random.AbstractRNG, parent::T, distr::MultivariatePolyCubeDistr) where { T <: MPolyRing }
  d = distr.deg
#  n = ngens(parent)
  vars = gens(parent)
  k = coefficient_ring(parent)
  poly = zero(parent)
  for i in 1:distr.nterms
    poly += prod([x^rand(RNG,0:d)  for x in vars]) * RAND(RNG, k, distr.CoeffDistr)
  end
  return poly
end

function DistrList(::T) where { T <: MPolyRing }
  return [MultivariatePolyCubeDistr];
end

# ------------------------------------------------------------------
# ----- (uniform, dense) Matrix -----

struct UnifDenseMatrixDistr <: AbstractRndDistr
  nr::Int  # number of rows
  nc::Int  # number of cols
  EntryDistr::AbstractRndDistr
  ApproxCardinality::Int64

  function UnifDenseMatrixDistr(r::Int, c::Int, entries::AbstractRndDistr)
    (r >= 0 && c >= 0) || error("Matrix dimensions must be non-negative")
    n = ApproxCard(entries)
    if n == 1
      card = 1
    elseif r*c >= 63
      card = typemax(Int64)
    else
      card = clamp(ZZ(n)^(r*c), 0, typemax(Int64))
      card = Int64(card)
    end
    return new(r,c, entries, card)
  end
end

function ContainsZero(D::UnifDenseMatrixDistr)
  return ContainsZero(D.EntryDistr)
end

function RAND(RNG::Random.AbstractRNG, parent::T, distr::UnifDenseMatrixDistr) where { T <: MatSpace }
  r = distr.nr
  c = distr.nc
  (nrows(parent) == r && ncols(parent) == c) || error("Incompatible dimensions")
  M = zero(parent)
  (r > 0 && c > 0) || return M
  for i in 1:r
    for j in 1:c
      M[i,j] = RAND(RNG, base_ring(parent), distr.EntryDistr)
    end
  end
  return M
end

# ------------------------------------------------------------------
# ----- (unimodular) Matrix -----

struct UnimodularMatrixDistr <: AbstractRndDistr
  n::Int  # number of rows/cols
  num_iters::Int
  EntryDistr::AbstractRndDistr
  ApproxCardinality::Int64

  function UnimodularMatrixDistr(n::Int, num_iters::Int, EntryDistr::AbstractRndDistr)
    (n >= 1) || error("Matrix dimensions must be positive")
    (num_iters >= 1) || error("Number of iteations must be positive")
    c = ApproxCard(EntryDistr)
    card = clamp(ZZ(n*c)^num_iters, 0, typemax(Int64))  # !!! UNDERESTIMATE !!!
    card = Int64(card)
    return new(n, num_iters, EntryDistr, card)
  end
end

function ContainsZero(D::UnifDenseMatrixDistr)
  return false;
end

function RAND(RNG::Random.AbstractRNG, parent::T, distr::UnimodularMatrixDistr) where { T <: MatSpace }
  n = distr.n
  num_iters = distr.num_iters
  M = parent(1)
  for k in 1:num_iters
    # Pick a pair of (distinct) row indexes: i, j
    i = rand(RNG, 1:n)
    j = rand(RNG, 1:n-1); if (j==i)  j = n; end
    add_row!(M, PostFilter(!is_zero, RNG, base_ring(parent), distr.EntryDistr), i, j)
  end
  return M
end

function DistrList(::T) where { T <: MatSpace }
  return [UnifDenseMatrixDistr, UnimodularMatrixDistr];
end

##################################################################

#function RAND(RNG::Random.AbstractRNG, parent, D::AbstractRndDistr)
#end

function RAND(parent, D::AbstractRndDistr)
  return RAND(Random.default_rng(), parent, D)
end


function RAND_(parent)
  return RAND(parent, DefaultDistribution(parent))
end

# ==================================================================

# Default distributions

function DefaultDistribution(::T) where {T <: FinField }
  return UniformDistr()
end

function DefaultDistribution(::T) where {T <: zzModRing }
  return UniformDistr()
end

function DefaultDistribution(::T) where {T <: ZZModRing }
  return UniformDistr()
end

function DefaultDistribution(::ZZRing)
  return UniformIntegerDistr(-3,3)
end

function DefaultDistribution(::UniformIntegerDistr)
  return UniformIntegerDistr(-3,3)
end

function DefaultDistribution(::QQField)
  return RationalDistr(DefaultDistribution(ZZ), UniformIntegerDistr(1,1))
end

function DefaultDistribution(::RationalDistr)
  return RationalDistr(DefaultDistribution(ZZ), UniformIntegerDistr(1,1))
end

function DefaultDistribution(P::T) where { T <: PolyRing }
  return UnivariatePolyDistr(2, DefaultDistribution(coefficient_ring(P)))
end

function DefaultDistribution(D::UnivariatePolyDistr)
  return UnivariatePolyDistr(2, DefaultDistribution(D.CoeffDistr))
end

function DefaultDistribution(P::T) where { T <: MPolyRing }
  return MultivariatePolyCubeDistr(2, 2, DefaultDistribution(coefficient_ring(P)))
end

function DefaultDistribution(D::MultivariatePolyCubeDistr)
  return MultivariatePolyCubeDistr(2, 2, DefaultDistribution(D.CoeffDistr))
end

function DefaultDistribution(S::T) where { T <: MatSpace }
  return UnifDenseMatrixDistr(nrows(S), ncols(S), DefaultDistribution(base_ring(S)))
end

function DefaultDistribution(D::UnifDenseMatrixDistr)
  return UnifDenseMatrixDistr(D.nr, D.nc, DefaultDistribution(D.EntryDistr))
end

# ------------------------------------------------------------------

# Widening

function WidenDistr(D::UniformDistr)
  return D  # do nothing
end

function WidenDistr(D::UniformIntegerDistr)
  width = D.upb-D.lwb
  if D.lwb < 0 && D.upb > 0
    return UniformIntegerDistr(2*D.lwb, 2*D.upb)
  end
  if D.lwb >= 0
    return UniformIntegerDistr(D.lwb, D.upb+width+1)
  else
    return UniformIntegerDistr(D.lwb-width-1, D.upb)
  end
end

function WidenDistr(D::RationalDistr)
  WidthNumer = width(D.NumerDistr)
  WidthDenom = width(D.DenomDistr)
  if 2*WidthDenom >= WidthNumer
    return RationalDistr(WidenDistr(D.NumerDistr), UniformIntegerDistr(1,1))
  end
  return RationalDistr(D.NumerDistr, WidenDistr(D.DenomDistr))
end

function WidenDistr(D::UnivariatePolyDistr)
  # Next line is boolean HEURISTIC whether to incr deg or CoeffDistr
  IncrDeg = (ApproxCard(DefaultDistribution(D.CoeffDistr)) < ApproxCard(D.CoeffDistr))
  if IncrDeg
    return UnivariatePolyDistr(1+D.deg, DefaultDistribution(D.CoeffDistr))
  end
  return UnivariatePolyDistr(D.deg, WidenDistr(D.CoeffDistr))
end

function WidenDistr(D::MultivariatePolyCubeDistr)
  # Next line is boolean HEURISTIC whether to incr deg or CoeffDistr
  if D.nterms <= D.deg^2
    return MultivariatePolyCubeDistr(D.deg, 1+D.nterms, DefaultDistribution(D.CoeffDistr))
  end
  IncrDeg = (ApproxCard(DefaultDistribution(D.CoeffDistr)) < ApproxCard(D.CoeffDistr))
  if IncrDeg
    return MultivariatePolyCubeDistr(1+D.deg, 2, DefaultDistribution(D.CoeffDistr))
  end
  return MultivariatePolyCubeDistr(D.deg, 2, WidenDistr(D.CoeffDistr))
end

function WidenDistr(D::UnifDenseMatrixDistr)
  return UnifDenseMatrixDistr(D.nr, D.nc, WidenDistr(D.EntryDistr))
end
