export LatticeWithIsometry, LWIType

@attributes mutable struct LatticeWithIsometry
  Lb::ZLat
  f::fmpq_mat
  f_ambient::fmpq_mat
  n::IntExt

  function LatticeWithIsometry(Lb::ZLat, f::fmpq_mat, f_ambient::fmpq_mat, n::IntExt)
    z = new()
    z.Lb = Lb
    z.f = f
    z.f_ambient = f_ambient
    z.n = n
    return z
  end
end

