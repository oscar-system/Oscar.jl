
########################################################################
# Rings of regular functions on Zariski open sets of affine schemes    #
########################################################################
@doc raw"""
    AffineSchemeOpenSubschemeRing{AffineSchemeType, OpenType}

The ring of regular functions ``𝒪(X, U)`` on an open subset ``U`` of an
affine scheme ``X``.

 * `AffineSchemeType` is the type of the affine scheme ``X`` on which
this sheaf is defined;
 * `OpenType` is the type of the (Zariski) open subsets of ``U``.
"""
mutable struct AffineSchemeOpenSubschemeRing{AffineSchemeType, OpenType} <: Ring
  scheme::AffineSchemeType
  domain::OpenType

  function AffineSchemeOpenSubschemeRing(
      X::AffineSchemeType,
      U::OpenType;
      check::Bool=true
    ) where {AffineSchemeType<:AbsAffineScheme, OpenType<:AffineSchemeOpenSubscheme}
    @check is_subscheme(U, X) "open set does not lay in the scheme"
    return new{AffineSchemeType, OpenType}(X, U)
  end
end

########################################################################
# Elements of AffineSchemeOpenSubschemeRings                                            #
########################################################################
@doc raw"""
    AffineSchemeOpenSubschemeRingElem{AffineSchemeOpenSubschemeType}

An element ``f ∈ 𝒪(X, U)`` of the ring of regular functions on
an open set ``U`` of an affine scheme ``X``.

The type parameter `AffineSchemeOpenSubschemeType` is the type of the open set
``U`` of ``X``.
"""
mutable struct AffineSchemeOpenSubschemeRingElem{
      AffineSchemeOpenSubschemeRingType<:AffineSchemeOpenSubschemeRing
    } <: RingElem
  parent::AffineSchemeOpenSubschemeRingType
  restrictions::Vector{<:RingElem}

  function AffineSchemeOpenSubschemeRingElem(
      R::AffineSchemeOpenSubschemeRingType,
      f::Vector{<:RingElem};
      check::Bool=true
    ) where {
        AffineSchemeOpenSubschemeRingType<:AffineSchemeOpenSubschemeRing
    }
    n = length(f)
    U = domain(R)
    n == length(affine_patches(U)) || error("the number of restrictions does not coincide with the number of affine patches")
    g = [OO(U[i])(f[i]) for i in 1:n] # will throw if conversion is not possible
    @check begin
      for i in 1:n-1
        for j in i+1:n
          W = U[i,j]
          OO(W)(f[i], check=false) == OO(W)(f[j], check=false) || error("elements are not compatible on overlap")
        end
      end
    end
    return new{AffineSchemeOpenSubschemeRingType}(R, g)
  end
end

