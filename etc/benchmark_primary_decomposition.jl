# Benchmark of the native primary decomposition algorithms in
# `src/Rings/groebner/primary-decomposition` against the corresponding procedures
# of Singular's `primdec.lib`. Run with
#   julia --project=. etc/benchmark_primary_decomposition.jl
using Oscar

# Best wall-clock time in seconds of `f` over `reps` runs after one warm-up run.
function best_time(f; reps::Int=3)
  f()
  return minimum(@elapsed(f()) for _ in 1:reps)
end

function cyclic(K, n)
  R, x = polynomial_ring(K, n, :x)
  gs = [sum(prod(x[mod1(i + j, n)] for j in 0:k-1) for i in 1:n) for k in 1:n-1]
  return ideal(R, push!(gs, prod(x) - 1))
end

function katsura(K, n)
  R, x = polynomial_ring(K, n + 1, :x)
  u(i) = abs(i) <= n ? x[abs(i)+1] : zero(R)
  gs = [sum(u(l) * u(m - l) for l in -n:n) - u(m) for m in 0:n-1]
  return ideal(R, push!(gs, sum(u(l) for l in -n:n) - 1))
end

function cases(K)
  R, (x, y, z, w) = polynomial_ring(K, [:x, :y, :z, :w])
  return [
    ("embedded (x^2,xy)", ideal(R, [x^2, x*y])),
    ("cone", ideal(R, [x^2*y, x*y^2, x*z^2, z^3*y])),
    ("embedded 2", ideal(R, [x^2 - y*z, y^3 - z^2*x, x*y*z])),
    ("sphere/plane", ideal(R, [(x^2+y^2+z^2-1)*(x-y), (x^2+y^2+z^2-1)*z])),
    ("algebraic", ideal(R, [(x^2-2)*(y^2+1), (x^2-2)*z, (y^2+1)*z^2])),
    ("lines", intersect(ideal(R, [x, y]), ideal(R, [y, z])^2, ideal(R, [x-z, y-1]))),
    ("points", ideal(R, [x*(x-1)*(x-2)*(x-3), y - x^2, z - x*y, w - x*z])),
    ("cyclic-3", cyclic(K, 3)),
    ("cyclic-4", cyclic(K, 4)),
    ("katsura-3", katsura(K, 3)),
  ]
end

function singular_call(f, I)
  S = Oscar.singular_polynomial_ring(I)
  G = Oscar.singular_generators(I)
  return () -> f(S, G)
end

const L = Oscar.Singular.LibPrimdec

for K in (QQ, GF(32003))
  println("\n== base field: ", K)
  println(rpad("example", 20), rpad("task", 12), rpad("native [s]", 14), rpad("Singular [s]", 14), "ratio")
  for (name, I) in cases(K)
    tasks = [
      ("minAss", () -> Oscar.minimal_associated_primes_gtz(I), singular_call(L.minAssGTZ, I)),
      ("primdecGTZ", () -> Oscar.primary_decomposition_gtz(I), singular_call(L.primdecGTZ, I)),
      ("primdecSY", () -> Oscar.primary_decomposition_sy(I), singular_call(L.primdecSY, I)),
      ("radical", () -> Oscar.radical_kl(I), singular_call(L.radical, I)),
      ("equidim", () -> Oscar.equidim(I), singular_call(L.equidim, I)),
    ]
    for (task, nat, sing) in tasks
      tn = try best_time(nat) catch e; NaN end
      ts = try best_time(sing) catch e; NaN end
      println(rpad(name, 20), rpad(task, 12), rpad(round(tn; sigdigits=3), 14),
              rpad(round(ts; sigdigits=3), 14), round(tn / ts; sigdigits=3))
    end
  end
end
