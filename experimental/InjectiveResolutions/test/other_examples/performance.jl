using Profile
using Statistics: median

# helper: run f() n times after one warm-up call, report times and allocations
function bench(label, f; samples=3)
  f()  # warm-up (compile + first-run overhead)
  times    = Float64[]
  gctimes  = Float64[]
  allocs   = Int[]
  bytes    = Int[]
  for _ in 1:samples
    # NOTE: deliberately NOT calling GC.gc() here. A forced GC before each
    # sample empties the heap and prevents collections from happening during
    # the timed run, making `stats.gctime` artificially 0. Letting the heap
    # accumulate gives a more realistic picture of GC overhead.
    stats = @timed f()
    push!(times,   stats.time)
    push!(gctimes, stats.gctime)
    push!(allocs,  Base.gc_alloc_count(stats.gcstats))
    push!(bytes,   stats.bytes)
  end
  gc_pct = 100 * median(gctimes) / median(times)
  println("\n=== $label ===")
  println("  samples : ", samples)
  println("  min     : ", round(minimum(times); digits=3), " s")
  println("  median  : ", round(median(times); digits=3), " s")
  println("  max     : ", round(maximum(times); digits=3), " s")
  println("  gc      : ", round(median(gctimes); digits=3), " s  (", round(gc_pct; digits=1), "%)")
  println("  alloc   : ", median(allocs), " allocations")
  println("  memory  : ", round(median(bytes) / 2^20; digits=1), " MiB")
end

# ---------------------------------------------------------------------------
# test 1: ZZ^2 cone, 2 generators
#   previous baseline (raw @time): 2.55 s
# ---------------------------------------------------------------------------
kQ = monoid_algebra([[1,0],[1,1]], QQ)
KQ = quotient_ring_as_module(ideal(kQ, []))
@time inj_res = injective_resolution(KQ, 3)
bench("test 1: monoid_algebra([[1,0],[1,1]]), i=3", () -> injective_resolution(KQ, 3); samples=5)

# === test 1: monoid_algebra([[1,0],[1,1]]), i=3 ===
#   samples : 5
#   min     : 1.103 s
#   median  : 1.173 s
#   max     : 1.324 s
#   alloc   : 610959.0 allocations
#   memory  : 18.5 MiB

# ---------------------------------------------------------------------------
# test 2: ZZ^2 cone, 3 generators
#   previous baseline (raw @time): 82.09 s
# ---------------------------------------------------------------------------
kQ = monoid_algebra([[1,0],[1,1],[1,2]], QQ)
KQ = quotient_ring_as_module(ideal(kQ, []))
@time inj_res = injective_resolution(KQ, 3)
bench("test 2: monoid_algebra([[1,0],[1,1],[1,2]]), i=3", () -> injective_resolution(KQ, 3); samples=2)

# === test 2: monoid_algebra([[1,0],[1,1],[1,2]]), i=3 ===
#   samples : 2
#   min     : 2.09 s
#   median  : 2.12 s
#   max     : 2.149 s
#   alloc   : 3.107115e6 allocations
#   memory  : 77.4 MiB

# ---------------------------------------------------------------------------
# test 3: ZZ^2 cone, 4 generators (slowest)
#   previous baseline (raw @time): 191.84 s
# ---------------------------------------------------------------------------
kQ = monoid_algebra([[1,0],[1,1],[1,2],[1,3]], QQ)
KQ = quotient_ring_as_module(ideal(kQ, []))

println("\n=== test 3: monoid_algebra([[1,0],[1,1],[1,2],[1,3]]), i=3 ===")
# only one sample — this is multi-minute
@time inj_res = injective_resolution(KQ, 3)
# 150.468644 seconds (264.47 M allocations: 6.836 GiB, 55.40% gc time)

# ---------------------------------------------------------------------------
# Profile pass on test 2 to locate hotspots
#
# Builds a flat report sorted by self-time, filtered to frames with at least
# `mincount` samples. The C=false flag hides C/runtime frames so the output
# focuses on Julia source lines.
#
# For an interactive flamegraph instead, install ProfileView:
#   using Pkg; Pkg.add("ProfileView")
#   using ProfileView; ProfileView.view()
# (or in VS Code, use the `@profview` macro from the Julia extension)
# ---------------------------------------------------------------------------
kQ = monoid_algebra([[1,0],[1,1],[1,2]], QQ)
KQ = quotient_ring_as_module(ideal(kQ, []))

println("\n=== profile pass on test 2 ===")
Profile.clear()
@profile injective_resolution(KQ, 3)

println("\n--- flat profile (top self-time hotspots) ---")
Profile.print(format=:flat, sortedby=:count, mincount=100, C=false)

println("\n--- tree profile (call structure, top 30 lines) ---")
Profile.print(format=:tree, maxdepth=20, mincount=200, C=false)
