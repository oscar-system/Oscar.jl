using Oscar
using GroebnerWalk
using BenchmarkTools

include("cyclic.jl")
include("katsura.jl")
include("agk.jl")
include("tran3.3.jl")
include("newellp1.jl")

function benchmark(
  io,
  name::String,
  I::Ideal,
  target::MonomialOrdering,
  start::MonomialOrdering
)
  print(io, name, ","); flush(io)
  # t = @belapsed groebner_walk($I, $target, $start; algorithm=:standard) seconds=20000 samples=5
  t = @elapsed groebner_walk(I, target, start; algorithm=:standard)
  print(io, t, ","); flush(io)
  # t = @belapsed groebner_walk($I, $target, $start; algorithm=:generic) seconds=20000 samples=5
  t = @elapsed groebner_walk(I, target, start; algorithm=:generic)
  print(io, t, ","); flush(io)
  # t = @belapsed groebner_basis($I; ordering=$target) seconds=20000 samples=5
  t = @elapsed groebner_basis(I; ordering=target)
  println(io, t); flush(io)
end

p = 11863279
Fp = GF(p)
open("results.csv", "a") do io
  R, (x, y) = QQ[:x, :y]
  I = ideal([y^4 + x^3 - x^2 + x, x^4])
  benchmark(io, "simple", I, lex(R), default_ordering(R))
  
  # benchmark(io, "cyclic5-QQ", cyclic5(QQ)...)
  # benchmark(io, "cyclic5-Fp", cyclic5(Fp)...)

  # benchmark(io, "cyclic6-QQ", cyclic6(QQ)...)
  # benchmark(io, "cyclic6-Fp", cyclic6(Fp)...)

  # benchmark(io, "cyclic7-QQ", cyclic7(QQ)...)
  # benchmark(io, "cyclic7-Fp", cyclic7(Fp)...)

  benchmark(io, "agk4-QQ", agk4(QQ)...)
  benchmark(io, "agk4-Fp", agk4(Fp)...)

  # benchmark(io, "katsura6-QQ", katsura6(QQ)...)
  # benchmark(io, "katsura6-Fp", katsura6(Fp)...)

  benchmark(io, "tran3.3-QQ", tran33(QQ)...)
  benchmark(io, "tran3.3-Fp", tran33(Fp)...)

  benchmark(io, "newellp1-QQ", newellp1(QQ)...)
  benchmark(io, "newellp1-Fp", newellp1(Fp)...)

  # R, (x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19, x20, x21, x22, x23, x24, x25, x26, x27, x28, x29, x30, x31, x32) = QQ[
  #   :x1, :x2, :x3, :x4, :x5, :x6, :x7, :x8, :x9, :x10, :x11, :x12, :x13, :x14, :x15, :x16, :x17, :x18, :x19, :x20, :x21, :x22, :x23, :x24, :x25, :x26, :x27, :x28, :x29, :x30, :x31, :x32
  # ]
  # F = [
  #   -x23*x32+x24*x31, -x22*x32+x24*x30, -x22*x31+x23*x30, -x21*x32+x24*x29, -x21*x31+x23*x29, -x21*x30+x22*x29, -x12*x32+x16*x28,  -x19*x28+x20*x27,
  #   -x11*x31+x15*x27, -x18*x28+x20*x26, -x18*x27+x19*x26, -x10*x30+x14*x26, -x17*x28+x20*x25, -x17*x27+x19*x25, -x17*x26+x18*x25, -x9*x29+x13*x25, x20*x8-x24*x4,
  #   -x17*x20-x17*x24-2*x17*x28-x17*x32+x18*x19+x18*x23+2*x18*x27+x18*x31+x19*x22+x19*x30-x20*x21-x20*x29-x21*x24-x21*x28-2*x21*x32+x22*x23+x22*x27+2*x22*x31+x23*x26-x24*x25-x25*x28-x25*x32+x26*x27+x26*x31+x27*x30-x28*x29-x29*x32+x30*x31,
  #   x19*x7-x23*x3, x18*x6-x2*x22, -x1*x21+x17*x5
  # ]
  # benchmark(io, "bayes148", ideal(R, F), lex(R), default_ordering(R))
end