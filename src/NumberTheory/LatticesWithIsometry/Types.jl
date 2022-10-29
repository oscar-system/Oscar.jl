export LatticeWithIsometry

@attributes mutable struct LatticeWithIsometry
  Lb::ZLat
  f::fmpq_mat
  f_ambient::fmpq_mat
  n::Integer

  function LatticeWithIsometry(Lb::ZLat, f::fmpq_mat, f_ambient::fmpq_mat, n::Integer)
    z = new()
    z.Lb = Lb
    z.f = f
    z.f_ambient = f_ambient
    z.n = n
    return z
  end
end

