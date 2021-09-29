mutable struct InvRing{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}
  field::FldT
  poly_ring::PolyRingT

  group::GrpT
  action::Vector{ActionT}
  action_singular::Vector{SingularActionT}

  modular::Bool

  primary::Vector{PolyElemT}
  secondary::Vector{PolyElemT}
  irreducible_secondary::Vector{PolyElemT}
  fundamental::Vector{PolyElemT}

  reynolds_operator::MapFromFunc{PolyRingT, PolyRingT}

  molien_series::Generic.Frac{fmpq_poly}

  # Cache some stuff on the Singular side
  # (possibly removed at some point)
  reynolds_singular::Singular.smatrix
  molien_singular::Singular.smatrix
  primary_singular # the type is different depending on the characteristic...

  function InvRing(K::FldT, G::GrpT, action::Vector{ActionT}) where {FldT <: Field, GrpT <: AbstractAlgebra.Group, ActionT}
    n = degree(G)
    R, = grade(PolynomialRing(K, "x" => 1:n, cached = false)[1], ones(Int, n))
    R_sing = singular_ring(R)
    action_singular = identity.([change_base_ring(R_sing, g) for g in action])
    PolyRingT = typeof(R)
    PolyElemT = elem_type(R)
    SingularActionT = eltype(action_singular)
    z = new{FldT, GrpT, PolyElemT, PolyRingT, ActionT, SingularActionT}()
    z.field = K
    z.poly_ring = R
    z.group = G
    z.action = action
    z.action_singular = action_singular
    z.modular = true
    if iszero(characteristic(K))
      z.modular = false
    else
      if !iszero(mod(order(G), characteristic(K)))
        z.modular = false
      end
    end
    return z
  end
end

struct AllMonomials{PolyRingT}
  R::PolyRingT
  d::Int

  function AllMonomials{PolyRingT}(R::PolyRingT, d::Int) where PolyRingT
    @assert d >= 0
    return new{PolyRingT}(R, d)
  end
end

struct InvRingBasisIterator{InvRingT, IteratorT, PolyElemT, MatrixT}
  R::InvRingT
  degree::Int
  dim::Int
  reynolds::Bool

  monomials::IteratorT
  monomials_collected::Vector{PolyElemT}
  kernel::MatrixT # used iff reynolds == false
end
