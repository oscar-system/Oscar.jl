module DiscLog

using Oscar
import Oscar.Nemo
import Oscar.Hecke

function __init__()
  Hecke.add_verbose_scope(:DiscLog)
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
# move to Hecke/Nemo/AA?
###############################################################
#TODO: consolidate with Hecke: is_primitive_root?

@attr Fac{fmpz} factored_order(K::FinField) = factor(size(K)-1)

function Oscar.is_primitive(a::FinFieldElem, f::Fac{fmpz} = factored_order(parent(a)))
  iszero(a) && return false
  n = size(parent(a))-1
  for p = keys(f.fac)
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

function disc_log(a::FinFieldElem, b::FinFieldElem)
  Nemo.check_parent(a, b)
  da = disc_log(a)
  db = disc_log(b)
  qm1 = size(parent(a))-1
  return (db*modinv(da, qm1)) % qm1
end

function disc_log(a::T) where {T <: FinFieldElem}
  K = parent(a)
  res = Vector{Tuple{fmpz, fmpz}}()
  qm1 = size(K)-1
  g = generator(K)
  for (p, k) = factored_order(K)
    r = divexact(qm1, p^k)
    push!(res, (fmpz(p)^k, Hecke.disc_log_ph(g^r, a^r, fmpz(p), k)))
  end
  return crt([x[2] for x = res], [x[1] for x= res])
end

export disc_log, generator

end # module DiscLog

export generator, disc_log

