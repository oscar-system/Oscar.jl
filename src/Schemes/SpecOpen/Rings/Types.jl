export SpecOpenRing
export SpecOpenRingElem

########################################################################
# Rings of regular functions on Zariski open sets of affine schemes    #
########################################################################
@Markdown.doc """
    SpecOpenRing{SpecType, OpenType}

The ring of regular functions ``𝒪(X, U)`` on an open subset ``U`` of an
affine scheme ``X``.

 * `SpecType` is the type of the affine scheme ``X`` on which
this sheaf is defined;
 * `OpenType` is the type of the (Zariski) open subsets of ``U``.
"""
mutable struct SpecOpenRing{SpecType, OpenType} <: Ring
  scheme::SpecType
  domain::OpenType

  function SpecOpenRing(
      X::SpecType,
      U::OpenType
    ) where {SpecType<:AbsSpec, OpenType<:SpecOpen}
    issubset(U, X) || error("open set does not lay in the scheme")
    return new{SpecType, OpenType}(X, U)
  end
end

########################################################################
# Elements of SpecOpenRings                                            #
########################################################################
@Markdown.doc """
    SpecOpenRingElem{SpecOpenType}

An element ``f ∈ 𝒪(X, U)`` of the ring of regular functions on
an open set ``U`` of an affine scheme ``X``.

The type parameter `SpecOpenType` is the type of the open set
``U`` of ``X``.
"""
mutable struct SpecOpenRingElem{
      SpecOpenRingType<:SpecOpenRing
    } <: RingElem
  parent::SpecOpenRingType
  restrictions::Vector{<:RingElem}

  function SpecOpenRingElem(
      R::SpecOpenRingType,
      f::Vector{<:RingElem};
      check::Bool=true
    ) where {
        SpecOpenRingType<:SpecOpenRing
    }
    n = length(f)
    U = domain(R)
    n == length(affine_patches(U)) || error("the number of restrictions does not coincide with the number of affine patches")
    g = [OO(U[i])(f[i]) for i in 1:n] # will throw if conversion is not possible
    if check
      for i in 1:n-1
        for j in i+1:n
          W = U[i,j]
          OO(W)(f[i], check=false) == OO(W)(f[j], check=false) || error("elements are not compatible on overlap")
        end
      end
    end
    return new{SpecOpenRingType}(R, g)
  end
end

