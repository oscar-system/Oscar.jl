########################################################################
# Printing                                                             #
########################################################################

# For the printing, we describe the domain and codomain with coordinates
# for the ambient space and the description of the complement (and we do not
# forget to avoid printing again the coordinates by using the `false`) argument
# in both the show function for `domain(f)` and `codomain(f)`).
function Base.show(io::IO, ::MIME"text/plain", f::AffineSchemeOpenSubschemeMor)
  io = pretty(io)
  X = domain(f)
  cX = ambient_coordinates(X)
  Y = codomain(f)
  cY = ambient_coordinates(Y)
  co_str = String[]
  str = "["*join(cX, ", ")*"]"
  kX = length(str)
  push!(co_str, str)
  str = "["*join(cY, ", ")*"]"
  kY = length(str)
  push!(co_str, str)
  k = max(length.(co_str)...)
  println(io, "Affine scheme open subscheme morphism")
  print(io, Indent(), "from ")
  print(io, co_str[1]*" "^(k-kX+2), Lowercase())
  show(IOContext(io, :show_coordinates => false), domain(f))
  println(io)
  print(io, "to   ", co_str[2]*" "^(k-kY+2), Lowercase())
  show(IOContext(io, :show_coordinates => false), codomain(f))
  mop = maps_on_patches(f)
  if length(mop) > 0
    println(io)
    print(io, Dedent(), "defined by the map")
    length(mop) > 1 && print(io, "s")
    print(io, Indent())
    for i in 1:length(mop)
      println(io, Lowercase())
      Base.show(io, MIME"text/plain"(), mop[i])
      if i != length(mop)
        println(io)
        print(io, "----------------------------------------------------------------------")
      end
    end
  end
  print(io, Dedent())
end

function Base.show(io::IO, f::AffineSchemeOpenSubschemeMor)
  if is_terse(io)
    print(io, "Affine scheme open subscheme morphism")
  else
    io = pretty(io)
    print(io, "Hom: ", Lowercase(), domain(f), " -> ", Lowercase(), codomain(f))
  end
end

########################################################################
# Composition of maps                                                  #
########################################################################
function compose(f::AffineSchemeOpenSubschemeMor, g::AffineSchemeOpenSubschemeMor; check::Bool=true)
  U = domain(f)
  Cf = codomain(f)
  V = domain(g)
  @check is_subscheme(Cf, V) "maps are not compatible"
  W = codomain(g)
  X = ambient_scheme(U)
  Y = ambient_scheme(V)
  Z = ambient_scheme(W)
  pb_coords = [pullback(f)(pullback(g)(OO(W)(x))) for x in gens(ambient_coordinate_ring(Z))]
  maps_on_patches = [morphism(A, Z, [restrict(h, A, check=check) for h in pb_coords], check=check) for A in affine_patches(U)]
  return AffineSchemeOpenSubschemeMor(U, W, maps_on_patches, check=check)
end

########################################################################
# Pullback of regular functions (deprecated)                           #
########################################################################
function pullback(f::AffineSchemeOpenSubschemeMor, a::RingElem)
  U = domain(f)
  X = ambient_scheme(U)
  V = codomain(f)
  Y = ambient_scheme(V)
  R = ambient_coordinate_ring(Y)
  parent(a) === R || error("element does not belong to the correct ring")
  pb_a = [pullback(f[i])(a) for i in 1:n_patches(U)]
  return AffineSchemeOpenSubschemeRingElem(AffineSchemeOpenSubschemeRing(X, U), pb_a)
end

########################################################################
# Equality test                                                        #
########################################################################
function ==(f::T, g::T) where {T<:AffineSchemeOpenSubschemeMor}
  (domain(f) == domain(g)) || return false
  (codomain(f) == codomain(g)) || return false
  Y = ambient_scheme(codomain(f))
  m = length(affine_patches(domain(f)))
  n = length(affine_patches(domain(g)))
  for i in 1:m
    for j in 1:n
      U = intersect(domain(f)[i], domain(g)[i])
      restrict(f[i], U, Y) == restrict(g[i], U, Y) || return false
    end
  end
  return true
end

########################################################################
# Preimages under AffineSchemeOpenSubschemeMors                                         #
########################################################################
function preimage(f::AffineSchemeOpenSubschemeMor, Z::AbsAffineScheme; check::Bool=true)
  U = domain(f)
  X = ambient_scheme(U)
  @check is_closed_embedding(Z, ambient_scheme(codomain(f))) "second argument must be closed in the codomain"
  n = length(affine_patches(U))
  pbZ = [preimage(f[i], Z) for i in 1:n]
  Y = X
  for K in pbZ
    Y = subscheme(Y, gens(modulus(underlying_quotient(OO(K)))))
  end
  return AffineSchemeOpenSubscheme(Y, [g for g in complement_equations(U) if !iszero(OO(Y)(g))])
end
function preimage(f::AffineSchemeOpenSubschemeMor, W::PrincipalOpenSubset; check::Bool=true)
  V = codomain(f)
  Y = ambient_scheme(V)
  Y === ambient_scheme(W) || error("second argument must be open in the ambient scheme of the domain of the morphism")
  h = complement_equation(W)
  pbh = pullback(f)(OO(codomain(f))(h))
  R = ambient_coordinate_ring(ambient_scheme(domain(f)))
  U = domain(f)
  X = ambient_scheme(U)
  I = ideal(R, one(R))
  for i in 1:n_patches(U)
    I = intersect(I, saturated_ideal(ideal(OO(U[i]), pbh[i])))
  end
  return intersect(U, AffineSchemeOpenSubscheme(X, I))
end


function preimage(f::AffineSchemeOpenSubschemeMor, V::AffineSchemeOpenSubscheme; check::Bool=true)
  U = domain(f)
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  I = ideal(R, one(R))
  for i in 1:n_patches(U)
    I = intersect(I, saturated_ideal(ideal(OO(U[i]), OO(U[i]).(pullback(f[i]).(complement_equations(V))))))
  end
  return intersect(U, AffineSchemeOpenSubscheme(X, I))
end

########################################################################
# Auxiliary methods                                                    #
########################################################################
function is_non_zero_divisor(f::RET, U::AffineSchemeOpenSubscheme) where {RET<:RingElem}
  return all(x->(is_non_zero_divisor(f, x)), affine_patches(U))
end

function find_non_zero_divisor(U::AffineSchemeOpenSubscheme)
  n = ngens(U)
  X = ambient_scheme(U)
  R = ambient_coordinate_ring(X)
  n == 0 && return zero(R)
  kk = base_ring(X)
  coeff = elem_type(kk)[rand(kk, 0:100) for i in 1:n]
  d = sum([coeff[i]*complement_equations(U)[i] for i in 1:n])
  while !is_non_zero_divisor(d, U)
    d = dot([rand(kk, 0:100) for i in 1:n], complement_equations(U))
  end
  return d
end

@doc raw"""
    generic_fractions(f::AffineSchemeOpenSubschemeMor)

Given a morphism ``f : U → V`` of Zariski open subsets ``U ⊂ X ⊂ 𝔸ᵐ`` and ``V ⊂ Y ⊂ 𝔸ⁿ``,
produce a tuple of fractions ``[a₁/b₁,…,aₙ/bₙ]`` such that ``f`` can be recovered
as the maximal extension of the rational map given by
```
   U ⊃ U' → 𝔸ⁿ,  x ↦ [a₁(x)/b₁(x),…,aₙ(x)/bₙ(x)]
```
where ``U'`` is the complement of the zero loci of the denominators ``bᵢ`` in ``U``.
In particular, this requires ``U'`` to be dense in ``U`` and this subset
is chosen at random.
"""
function generic_fractions(f::AffineSchemeOpenSubschemeMor)
  U = domain(f)
  X = ambient_scheme(U)
  V = codomain(f)
  Y = ambient_scheme(V)
  d = find_non_zero_divisor(U)
  W = hypersurface_complement(X, d)
  result = fraction.([restrict(pullback(f)(OO(V)(y)), W) for y in gens(ambient_coordinate_ring(Y))])
  return result
end

