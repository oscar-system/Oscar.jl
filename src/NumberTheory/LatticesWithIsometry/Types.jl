@attributes mutable struct LatticeWithIsometry{S, T} <: AbsLat{S}
  Lb::AbsLat{S}
  f::T
  GL
  Lh::HermLat
  n::Int

  function LatticeWithIsometry{S, T}(Lb::AbsLat{S}, f::T, GL, Lh::HermLat, n::Int) where {S, T}
    z = new{S, T}(Lb, f, GL, Lh, n)
    return z
  end

  function LatticeWithIsometry{S, T}(Lb::AbsLat{S}, f::T, n::Int) where {S, T}
    z = new{S, T}()
    z.Lb = Lb
    z.f = f
    z.n = n
    Lh = Hecke.inverse_trace_lattice(Lb, f)
    z.Lh = Lh
    return z
  end
end

lattice(Lf::LatticeWithIsometry) = Lf.Lb

isometry(Lf::LatticeWithIsometry) = Lf.f

hermitian_structure(Lf::LatticeWithIsometry) =  Lf.Lh

rank(Lf::LatticeWithIsometry) = rank(Lf.Lb)

degree(Lf::LatticeWithIsometry) = degree(Lf.Lb)

genus(Lf::LatticeWithIsometry) = genus(Lf.Lb)

minpoly(Lf::LatticeWithIsometry) = minpoly(isometry(Lf))

charpoly(Lf::LatticeWithIsometry) = charpoly(isometry(Lf))

image_centralizer_in_Oq(Lf::LatticeWithIsometry) = Lf.GL

order_isometry(Lf::LatticeWithIsometry) = Lf.n

