using Oscar

include("cyclic.jl")
include("katsura.jl")
include("agk.jl")
include("tran3.3.jl")

function benchmark(
  io,
  name::String,
  I::Ideal,
  target::MonomialOrdering,
  start::MonomialOrdering
)
  print(io, name, ","); flush(io)
  t = @elapsed groebner_walk($I, $target, $start; algorithm=:standard) 
  print(io, t, ","); flush(io)
  t = @elapsed groebner_walk($I, $target, $start; algorithm=:generic)
  print(io, t, ","); flush(io)
  t = @elapsed groebner_basis($I; ordering=$target)
  println(io, t); flush(io)
end

function print_header(io)
  print(io, "name,standard_walk,generic_walk,buchberger\n")
end

p = 11863279
Fp = GF(p)
open("results.csv", "a") do io
  R, (x, y) = QQ[:x, :y]
  I = ideal([y^4 + x^3 - x^2 + x, x^4])
  benchmark(io, "simple", I, lex(R), default_ordering(R))
  
  benchmark(io, "cyclic5-QQ", cyclic5(QQ)...)
  benchmark(io, "cyclic5-Fp", cyclic5(Fp)...)

  benchmark(io, "cyclic6-QQ", cyclic6(QQ)...)
  benchmark(io, "cyclic6-Fp", cyclic6(Fp)...)
  
  I = katsura(6)
  R = base_ring(I)
  benchmark(io, "katsura6-QQ", I, lex(R), default_ordering(R))

  Ip = map(gens(I)) do f
    change_coefficient_ring(Fp, f)
  end |> ideal
  Rp = base_ring(Ip)
  benchmark(io, "katsura6-Fp", Ip, lex(Rp), default_ordering(Rp))

  benchmark(io, "cyclic7-QQ", cyclic7(QQ)...)
  benchmark(io, "cyclic7-Fp", cyclic7(Fp)...)

  benchmark(io, "agk4-QQ", agk4(QQ)...)
  benchmark(io, "agk4-Fp", agk4(Fp)...)

  benchmark(io, "tran3.3-QQ", tran33(QQ)...)
  benchmark(io, "tran3.3-Fp", tran33(Fp)...)

  benchmark(io, "newellp1-QQ", Oscar.newell_patch_with_orderings(QQ)...)
  benchmark(io, "newellp1-Fp", Oscar.newell_patch_with_orderings(Fp)...)
end

