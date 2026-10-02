# Prototype date: 2026-04-14

# Ad hoc documentation (just as a comment currently)

# (BUGGY) PROTOTYPE for generating random values

# Summary of features for this prototype (in no special order)
# (1) explicit distribution objects
# (2) RAND1(parent, distrib)  or  RAND1(RNG, parent, distrib)
# (3) iterator for random values: use `first(iter)` to get values
# (4) NOT YET IMPLEMENTED: DefaultDistribution(parent)
# (5) NOT YET IMPLEMENTED: SmallDistr and WidenDistr(D)
# (6) NOT YET IMPLEMENTED: ApproxPopulationSize

# Function/struct and other names are _ad hoc_; similarly for fn prototypes.


# ==================================================================
# Some basic "distribution" object types

# "abstract nonsense" not needed?

abstract type AbstractRndDistr end

function ContainsZero(D::AbstractRndDistr) end

# function ApproxCard(D::AbstractRndDistr)
#   return D.ApproxCardinality
# end

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


# ------------------------------------------------------------------
# Iterator interface  (preliminary)
# Uncertain what the second return value should be in this case

abstract type AbstractRndIter end


function RAND(iter::AbstractRndIter)
end



# ------------------------------------------------------------------
# ----- Finite Fields -----
# To achieve "type stability", I need to handle the following cases separately
# fpField, fqPolyRepField, FpField, FqField/FqPolyRepField

abstract type FiniteFieldDistr <: AbstractRndDistr end  # not sure if this is a good idea -- impedes "discoverability"

struct FiniteFieldDistr_SmallPrimeSubfield
  p::Int  # prime

  function FiniteFieldDistr_SmallPrimeSubfield(p::Int)   # arg is characteristic of the field
    (is_prime(p)) || error("arg must be prime");
    return new(p)
  end
end

function ContainsZero(D::FiniteFieldDistr_SmallPrimeSubfield)
  return true
end

struct FiniteField_SmallPrimeSubfieldIter  <: AbstractRndIter
  parent::Union{fpField, fqPolyRepField}
  distr::FiniteFieldDistr_SmallPrimeSubfield
  rng::Random.AbstractRNG

  function FiniteField_SmallPrimeSubfieldIter(k::Union{fpField, fqPolyRepField}, rng::AbstractRNG=Random.default_rng())
    return new(k, FiniteFieldDistr_SmallPrimeSubfield(Int(characteristic(k))), rng);
  end
end


function RAND(iter::FiniteField_SmallPrimeSubfieldIter)
  p = iter.distr.p
  return iter.parent(rand(iter.rng, 0:p-1));
end

function RAND1(parent::Union{fpField,fqPolyRepField}, D::FiniteFieldDistr_SmallPrimeSubfield)
  p = size(parent);
  (p == D.p) || error("Incompatible characteristics: field and distribution")
  return parent(rand(0:p-1));
end
  
function RAND1(RNG::Random.AbstractRNG, parent::fpField, D::FiniteFieldDistr_SmallPrimeSubfield)
  p = size(parent);
  (p == D.p) || error("Incompatible characteristics: field and distribution")
  return parent(rand(rng, 0:p-1));
end


#---------------------------------

struct FiniteFieldDistr_LargePrimeSubfield
  p::ZZRingElem  # prime

  function FiniteFieldDistr_LargePrimeSubfield(p::ZZRingElem) # arg is characteristic of the field
    (is_probable_prime(p)) || error("arg must be prime");
    return new(p)
  end
end

function ContainsZero(D::FiniteFieldDistr_LargePrimeSubfield)
  return true
end

struct FiniteField_LargePrimeSubfieldIter  <: AbstractRndIter
  parent::Union{FpField, FqField}
  distr::FiniteFieldDistr_LargePrimeSubfield
  rng::Random.AbstractRNG

  function FiniteField_LargePrimeSubfieldIter(k::Union{FpField, FqField}, rng::AbstractRNG=Random.default_rng())
    return new(k, FiniteFieldDistr_LargePrimeSubfield(characteristic(k)), rng);
  end
end


function RAND(iter::FiniteField_LargePrimeSubfieldIter)
  p = iter.distr.p;
  return iter.parent(rand(iter.rng, 0:p-1));
end

function RAND1(parent::Union{FpField, FqField}, D::FiniteFieldDistr_LargePrimeSubfield)
  p = characteristic(parent);
  (p == D.p) || error("Incompatible characteristics: field and distribution")
  return parent(rand(0:p-1));
end
  
function RAND1(RNG::Random.AbstractRNG, parent::Union{FpField,FqField}, D::FiniteFieldDistr_LargePrimeSubfield)
  p = size(parent);
  (p == D.characteristic) || error("Incompatible characteristics: field and distribution")
  return parent(rand(rng, 0:p-1));
end


# -----------------------------------------------------------------------------

struct FiniteFieldDistr_Uniform
#  characteristic::Int

  function FiniteFieldDistr_Uniform(#=q::Int=#)  
#    (q > 1 && is_prime_power_with_data(p)[1]) || error("cardinality must be a prime power");
    return new(#=q=#)
  end
end

function ContainsZero(D::FiniteFieldDistr_Uniform)
  return true
end

struct FiniteField_UniformIter  <: AbstractRndIter
  parent::Union{fpField, fqPolyRepField, FpField, FqField}
  distr::FiniteFieldDistr_Uniform
  rng::Random.AbstractRNG

  function FiniteField_UniformIter(k::Union{fpField, fqPolyRepField, FpField, FqField}, D::FiniteFieldDistr_Uniform, rng::AbstractRNG=Random.default_rng())
    return new(k, FiniteFieldDistr_Uniform(#=Int(characteristic(k))=#), rng);
  end
end


function RAND(iter::FiniteField_UniformIter)
#  p = iter.distr.characteristic;
  return rand(iter.rng, iter.parent);
end

function RAND1(parent::Union{fpField,fqPolyRepField, FpField, FqField}, D::FiniteFieldDistr_Uniform)
#  p = size(parent);
#  (p == D.characteristic) || error("Incompatible characteristics: field and distribution")
  return rand(parent);
end
  
function RAND1(RNG::Random.AbstractRNG, parent::Union{fpField, fqPolyRepField, FpField, FqField}, D::FiniteFieldDistr_Uniform)
#  p = size(parent);
#  (p == D.characteristic) || error("Incompatible characteristics: field and distribution")
  return rand(rng, parent);
end


# #---------------------------------

# struct FiniteFieldDistr_LargeUniform
#   characteristic::ZZRingElem

#   function FiniteFieldDistr_LargeUniform(p::ZZRingElem)
#     (is_prime_power_with_data(p)[1]) || error("cardinality must be a prime power");
#     return new(p)
#   end
# end

# function ContainsZero(D::FiniteFieldDistr_LargeUniform)
#   return true
# end

# struct FiniteField_LargeUniformIter  <: AbstractRndIter
#   parent::Union{FpField, FqField}
#   distr::FiniteFieldDistr_LargeUniform
#   rng::Random.AbstractRNG

#   function FiniteField_LargeUniformIter(k::Union{FpField, FqField}, rng::AbstractRNG=Random.default_rng())
#     return new(k, FiniteFieldDistr_LargeUniform(characteristic(k)), rng);
#   end
# end


# function RAND(iter::FiniteField_LargeUniformIter)
#   p = iter.distr.characteristic;
#   return rand(iter.rng, iter.parent);
# end

# function RAND1(parent::Union{FpField, FqField}, D::FiniteFieldDistr_LargeUniform)
#   p = size(parent);
#   (p == D.characteristic) || error("Incompatible characteristics: field and distribution")
#   return rand(parent);
# end
  
# function RAND1(RNG::Random.AbstractRNG, parent::Union{FpField,FqField}, D::FiniteFieldDistr_LargeUniform)
#   p = size(parent);
#   (p == D.characteristic) || error("Incompatible characteristics: field and distribution")
#   return rand(rng, parent);
# end



# ==================================================================
# Integers

struct IntegerDistr_UniformRange
  range::ZZRingElemUnitRange
#  lo::ZZRingElem
#  hi::ZZRingElem

  function IntegerDistr_UniformRange(range::ZZRingElemUnitRange)
    is_empty(range) && error("Sampling range must be non-empty");
    return new(range);
  end
end

# convenience ctors
IntegerDistr_UniformRange(lwb::Nemo.IntegerUnion, upb::Nemo.IntegerUnion) = IntegerDistr_UniformRange(ZZ(lwb):ZZ(upb));
IntegerDistr_UniformRange(r::UnitRange) = IntegerDistr_UniformRange(ZZ(first(r)):ZZ(last(r)));



function ContainsZero(D::IntegerDistr_UniformRange)
  return (first(D.range) <= 0 && last(D.range) >= 0);
end

struct IntegerDistr_UniformRangeIter  <: AbstractRndIter
  parent::ZZRing
  distr::IntegerDistr_UniformRange
  rng::Random.AbstractRNG

  function IntegerDistr_UniformRangeIter(Z::ZZRing, D::IntegerDistr_UniformRange, rng::AbstractRNG=Random.default_rng())
    return new(Z, D, rng);
  end
end

function RAND(iter::IntegerDistr_UniformRangeIter)
  return iter.parent(rand(iter.rng, iter.distr.range));
end

function RAND1(parent::ZZRing, D::IntegerDistr_UniformRange)
  return parent(rand(D.range));
end
  
function RAND1(RNG::Random.AbstractRNG, parent::ZZRing, D::IntegerDistr_UniformRange)
  return parent(rand(rng, D.range));
end


# ==================================================================
# Rationals

struct RationalDistr_NumDenRange
  NumRange::ZZRingElemUnitRange;
  DenRange::ZZRingElemUnitRange;

  function RationalDistr_NumDenRange(; NumRange::ZZRingElemUnitRange, DenRange::ZZRingElemUnitRange)
    (is_empty(NumRange) || is_empty(DenRange)) && error("Sampling range must be non-empty");
    (first(DenRange) >= 1) || error("Denominator range must be positive");
    return new(NumRange, DenRange);
  end
end

# # convenience ctors
RationalDistr_NumDenRange(Num_lo::Nemo.IntegerUnion, Num_hi::Nemo.IntegerUnion, Den_lo::Nemo.IntegerUnion, Den_hi::Nemo.IntegerUnion) = RationalDistr_NumDenRange(; NumRange=ZZ(Num_lo):ZZ(Num_hi), DenRange=ZZ(Den_lo):ZZ(Den_hi));
RationalDistr_NumDenRange(NumRange::UnitRange, DenRange::UnitRange) = RationalDistr_NumDenRange(; NumRange=ZZ(first(NumRange)):ZZ(last(NumRange)), DenRange=ZZ(first(DenRange)):ZZ(last(DenRange)));


Base.show(io::IO, D::RationalDistr_NumDenRange) = print(io, "RationalDistr(NumRange=$(D.NumRange), DenRange=$(D.DenRange))");


function ContainsZero(D::RationalDistr_NumDenRange)
  return (first(D.NumRange) <= 0 && last(D.NumRange) >= 0);
end

struct RationalDistr_NumDenRangeIter  <: AbstractRndIter
  parent::QQField
  distr::RationalDistr_NumDenRange
  rng::Random.AbstractRNG

  function RationalDistr_NumDenRangeIter(Q::QQField, D::RationalDistr_NumDenRange, rng::AbstractRNG=Random.default_rng())
    return new(Q, D, rng);
  end
end

function RAND(iter::RationalDistr_NumDenRangeIter)
  n = rand(iter.rng, iter.distr.NumRange);
  d = rand(iter.rng, iter.distr.DenRange);
  return iter.parent(n//d);
end

function RAND1(parent::QQField, D::RationalDistr_NumDenRangeIter)
  n = rand(D.NumRange);
  d = rand(D.DenRange);
  return parent(n//d);
end
  
function RAND1(RNG::Random.AbstractRNG, parent::QQField, D::RationalDistr_NumDenRangeIter)
  n = rand(RNG, D.NumRange);
  d = rand(RNG, D.DenRange);
  return parent(n//d);
end


# ==================================================================
# Univariate polynomials

struct PolyDistr_Dense
  deg::Int64
  CoeffGenerator::AbstractRndIter

  function PolyDistr_Dense(; deg::Int64, CoeffGenerator::AbstractRndIter)
    (deg >= 0) || error("Degree must be non-negative");
    return new(deg, CoeffGenerator);
  end
end



function ContainsZero(D::PolyDistr_Dense)
  return ContainsZero(D.CoeffGenerator.distr);
end

struct PolyDistr_DenseIter  <: AbstractRndIter
  parent::PolyRing
  distr::PolyDistr_Dense
  rng::Random.AbstractRNG

  function PolyDistr_DenseIter(P::PolyRing, DistrParams::PolyDistr_Dense, rng::AbstractRNG=Random.default_rng())
    return new(P, DistrParams, rng);
  end
end

function GENERATE_SAMPLE(rng::AbstractRNG, parent::PolyRing, distr::PolyDistr_Dense)
  P = parent;
  d = distr.deg;
  x = gen(P);
  f = zero(P);
  for i in 0:d
    f += RAND(distr.CoeffGenerator)*x^i;
  end
  return f;
end


function RAND(iter::PolyDistr_DenseIter)
  return GENERATE_SAMPLE(iter.rng, iter.parent, iter.distr);
end

function RAND1(parent::PolyRing, D::PolyDistr_Dense)
  return GENERATE_SAMPLE(Random.default_rng(), parent, D);
end
  
function RAND1(RNG::Random.AbstractRNG, parent::PolyRing, D::PolyDistr_Dense)
  return GENERATE_SAMPLE(RNG, parent, D);
end



# ==================================================================
# Multivariate polynomials

struct MPolyDistr_Simplex
  deg::Int64
  NumTerms::Int64
  CoeffGenerator::AbstractRndIter

  function MPolyDistr_Simplex(; deg::Int64, NumTerms::Int64, CoeffGenerator::AbstractRndIter)
    (deg >= 0) || error("Degree must be non-negative");
    (NumTerms >= 0) || error("Number of terms must be non-negative");
    return new(deg, NumTerms, CoeffGenerator);
  end
end


function ContainsZero(D::MPolyDistr_Simplex)
  # not quite correct!
  return (D.NumTerms == 0 || ContainsZero(D.CoeffGenerator.distr));
end

struct MPolyDistr_SimplexIter  <: AbstractRndIter
  parent::MPolyRing
  distr::MPolyDistr_Simplex
  rng::Random.AbstractRNG

  function MPolyDistr_SimplexIter(P::MPolyRing, DistrParams::MPolyDistr_Simplex, rng::AbstractRNG=Random.default_rng())
    (base_ring(P) == DistrParams.CoeffGenerator.parent) || error("Parent mismatch");
    return new(P, DistrParams, rng);
  end
end

function GENERATE_SAMPLE(rng::AbstractRNG, parent::MPolyRing, distr::MPolyDistr_Simplex)
  is_zero(distr.NumTerms) && return zero(P);
  P = parent;
  n = ngens(P);
  d = distr.deg;
  MaxNumTerms = binomial(n+d-1,d);
  NumIters = min(distr.NumTerms, MaxNumTerms);
  x = gens(P);
  f = zero(P);
  for _ in 0:NumIters
    ordered_subset = sort(first(Random.randperm(n+d),n))
    exps = [0  for _ in 1:n];
    exps[1] = ordered_subset[1]-1;
    for k in 2:n
      exps[k] = ordered_subset[k]-ordered_subset[k-1]-1;
    end
    t = prod(x.^exps);   # Is there a better way?
    f += RAND(distr.CoeffGenerator)*t;
  end
  return f;
end

function RAND(iter::MPolyDistr_SimplexIter)
  return GENERATE_SAMPLE(iter.rng, iter.parent, iter.distr);
end

function RAND1(parent::MPolyRing, D::MPolyDistr_Simplex)
  return GENERATE_SAMPLE(Random.default_rng(), parent, distr);
end
  
function RAND1(RNG::Random.AbstractRNG, parent::MPolyRing, D::MPolyDistr_Simplex)
  return GENERATE_SAMPLE(RNG, parent, distr);
end


# =======================================================
# (Dense) Matrices

struct MatrixDistr_DenseUniform
   nr::Int  # number of rows
   nc::Int  # number of cols
   EntryGenerator::AbstractRndIter

  function MatrixDistr_DenseUniform(; nr::Int64, nc::Int64, EntryGenerator::AbstractRndIter)
    (nr >= 0 && nc >= 0) || error("Degree must be non-negative");
    return new(nr, nc, EntryGenerator);
  end
end



function ContainsZero(D::MatrixDistr_DenseUniform)
  # not quite correct!
  return (D.nr == 0 || D.nc == 0 || ContainsZero(D.EntryGenerator.distr));
end

struct MatrixDistr_DenseUniformIter  <: AbstractRndIter
  parent::MatSpace
  distr::MatrixDistr_DenseUniform
  rng::Random.AbstractRNG

  function MatrixDistr_DenseUniformIter(S::MatSpace, DistrParams::MatrixDistr_DenseUniform, rng::AbstractRNG=Random.default_rng())
    (base_ring(S) == DistrParams.EntryGenerator.parent) || error("Parent mismatch");
    (nrows(S) == DistrParams.nr && ncols(S) == DistrParams.nc) || error("Matrix dimension mismatch");
    return new(S, DistrParams, rng);
  end
end

function GENERATE_SAMPLE(rng::AbstractRNG, parent::MatSpace, distr::MatrixDistr_DenseUniform)
  (is_zero(distr.nr) || is_zero(distr.nc)) && return zero(parent);
  M = zero(parent);
  for i in 1:distr.nr
    for j in 1:distr.nc
      M[i,j] = RAND(distr.EntryGenerator);
    end;
  end;
  return M;
end

function RAND(iter::MatrixDistr_DenseUniformIter)
  return GENERATE_SAMPLE(iter.rng, iter.parent, iter.distr);
end

function RAND1(parent::MatrixSpace, D::MatrixDistr_DenseUniform)
  return GENERATE_SAMPLE(Random.default_rng(), parent, distr);
end
  
function RAND1(RNG::Random.AbstractRNG, parent::MatrixSpace, D::MatrixDistr_DenseUniform)
  return GENERATE_SAMPLE(RNG, parent, distr);
end




#==================================================================
2026-04-08  API musings

RAND_ITER(parent, distr)
RAND_ITER(parent, distr, RandomnessSource)

RAND1(parent)                   -- uses default distr
RAND1(parent, RandomnessSource) -- uses default distr

A simple implementation of RAND1 could create an iterator,
then use it to generate one randon value, then discard the
iterator.  Not efficient, but simple, and guarantees
"compatibility" with whatever the iterator does.

==================================================================#

#==================================================================
NOTE ABOUT widening/expanding

To explain what I mean by widening/expanding, consider the following
(pseudo-)code snippet for searching for a small random value satisfying
some nice property NiceProp.  The idea is to start with sampling from
a small population, and if that produced no witness, expand the
population and restart the search.

  distr = InitialSmallDistribution(...);
  NumTrials = 100;  # value depends on likelihood of NiceProp being satisfied
  found = false;
  while true
    for i in 1:NumTrials
      x = RAND(distr);
      if NiceProp(x)
        found = true;
        witness = x;
        break;
      end
    end
    if found  break; end
    distr = WidenDistr(distr);
end
 
The API needs to offer:
(1)  InitialSmallDistribution
     [distr over some small initial population]
(2)  WidenDistr or ExpandDistr
     [produce a new distr over a larger population]
(3)  IsExpandable
     [say whether a larger population exists, e.g. in a finite ring]

==================================================================#

#==================================================================
NOTE ABOUT "default" distr vs. "small" distr

The "default" distr is the one which is used when a "lazy" caller
does not want to create a distr object prior to asking for some
random samples.  For instance  RAND_DEFAULT(ZZRing) might produce
values in the range -99:99 (or maybe even just -9:9)

Let K be a non-prime finite field.  I believe we want the following
distinct behaviours:
(1)  default distr is uniform over the whole of K
(2)  small distr is uniform over the prime subfield of K
     [or possibly even a subset of the prime subfield, e.g. {-1,0,1}]

==================================================================#


#==================================================================
NOTES ABOUT discoverability

It is not quite clear to me who would use this feature and when.
In principle, the information would be easy to find in well-organized
documentation (perhaps searchable via some LLM interface?).

In the section about "API musings" there were 2 functions:
RAND_ITER and RAND1.  One could use Julia's TAB completion with
"RAND_ITER(parent" to find out the names of types of distr which are
compatible with the given parent.  This might require a little care
when defining RAND_ITER...

==================================================================#


#==================================================================
NOTES ABOUT estimating populatiom size

Claus Fieker expressed interest in being able to obtain an estimate
of the population size for the random samples.  Having such an
estimate allows the caller to write code which avoids processing
too many duplicate random values (without having to keep a list
of the values already processed).

A significant problem is that it may be difficult to obtain a good
estimate of the population size.  Already this is tricky for
rational numbers with numerators in a specified range, and denominators
in another specified (positive) range.  However see file
   2026-03-30-CountDistinctRationals.cocoa5
for some code which produces good lower and upper bounds for the
population size -- the code is tolerably fast provided one of the
ranges is not too big.  Both ranges must be of the form 1:UPB.
There is an asymptotic formula for the population size, but I'm
not sure how accurate it is for small ranges...

Here are some API/UI suggestions:
(1)  is_population_large
     returns boolean saying whether population size > 1000000
(2)  is_population_probably_large:
     like (1) but based on probably good guess of pop size
(3)  is_population_larger_than(SIZE)
     simlar to (1) but caller can specify the critical size
(4)  is_population_probably_larger_than(SIZE)
     like (3) but based on probably good guess of pop size
(5)  population_size_range
     return a range LWB:UPB such that population size is within the range
(6)  like (5) but based on probably good guess of pop size
(7)  ???

One problem is that for some distributions the best answer is
quite possibly "Sorry, I haven't got a clue!"  How does API handle this?

WHAT ARGUMENT(S) must be passed to the functions which estimate
the population size?  It is possible that the distr alone does not
contain enough information (e.g. if there is a distr of type
UniformOnPrimeSubfield)

==================================================================#
