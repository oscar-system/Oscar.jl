module DiscLog  # TODO: move to Hecke

using Oscar
import Oscar.Nemo
import Oscar.Hecke
import Oscar: disc_log, is_primitive

function __init__()
  Hecke.add_verbosity_scope(:DiscLog)
end

###############################################################
#
# disc_log for finite fields
# TODO
# use the Conway property if applicable
# store the dlog data on the field to re-use (DONE)
# do more intelligent (sub-exponential) algorithms
# (maybe not all, but at least some?)
# see Hecke/src/Misc/UnitsModM: for BS-GS, cache the BS Array
# sort the types...
# TODO: move to Hecke/Nemo/AA?
###############################################################
#TODO: consolidate with Hecke: is_primitive_root?

@attr Fac{ZZRingElem} factored_order(K::FinField) = factor(size(K)-1)

function is_primitive(a::FinFieldElem, f::Fac{ZZRingElem} = factored_order(parent(a)))
  iszero(a) && return false
  n = size(parent(a))-1
  for (p, _) in f
    if a^divexact(n, p) == 1
      return false
    end
  end
  return true
end

@attr elem_type(T) function generator(K::T) where {T <: FinField}
  #if isconway use gen(K)
  a = rand(K)
  while !is_primitive(a)
    a = rand(K)
  end
  return a
end

function element_of_given_order(K::T, o::ZZRingElem) where {T <: FinField}
  @req is_positive(o) "order must be positive"
  q = order(K)
  e = divexact(q-1, o)
  lp = [divexact(o, p) for (p,k) in factor(o)]
  while true
    a = rand(K)^e
    (!is_zero(a) && all(!isone(a^x) for x in lp)) && return a
  end
end

function disc_log(a::FinFieldElem, b::FinFieldElem)
  Nemo.check_parent(a, b)
  da = disc_log(a)
  db = disc_log(b)
  qm1 = size(parent(a))-1
  return (db*invmod(da, qm1)) % qm1
end

function disc_log(a::T) where {T <: FinFieldElem}
  K = parent(a)
  res = Vector{Tuple{ZZRingElem, ZZRingElem}}()
  qm1 = size(K)-1
  g = generator(K)
  for (p, k) = factored_order(K)
    r = divexact(qm1, p^k)
    push!(res, (ZZRingElem(p)^k, Hecke.disc_log_ph(g^r, a^r, ZZRingElem(p), k)))
  end
  return crt([x[2] for x = res], [x[1] for x= res])
end

end # module DiscLog
