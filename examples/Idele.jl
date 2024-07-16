module Idele

using Oscar

mutable struct IdeleParent
  k::AbsSimpleNumField
  mG::Map # AutGrp -> Automorohisms
  S::Vector{AbsNumFieldOrderIdeal} # for each prime number ONE ideal above
  C::Vector{Map} # the completions at S
  L::Vector{Map} # the mult. group map at C

  #for P in S the modules used actually is
  #    Ind_G_p^G L[P]
  #        = sum L[P] otimes s_i
  # (for s_i a fixed system of coset reps G//G_P)
  # L[P] otimes s_i "should be" the completion data at P^s_i - one of the other ideals
  # should be L[P] ni l -> C[P] -> k -> inv(s_i)(..) to get a proper rep in k
  # completion at P^s is C[P] but with the map twisted by s

  mU::Map #S-unit group map
  M::FinGenAbGroup  # the big module

  function IdeleParent()
    return new()
  end
end

mutable struct Idele
  parent::IdeleParent
  m::FinGenAbGroupElem #in parent.M
end

function support(a::Idele)
  #the full galois orbit of parent.S
end

function getindex(a::Idele, P::AbsNumFieldOrderIdeal)
  #element at place P as an element in the completion
  #needs to find the "correct" P in parent.S, the (index) of the coset
  #and return the component of the induced module
end

function getindex(a::Idele, p) #real embedding
  # the unit representing this: find the correct component...
end

function getindex(a::Idele, p) #complex embedding
  # the unit representing this + a vector with the lambda exponents
end

function disc_exp(a::Idele)
  # a dict mapping primes/places to s.th. in the number field?
end

function disc_log()
end


end #module
