orbit_counts = MSet{Int}();

for k in 1:92
  J = johnson_solid(k)
  Aut = automorphism_group(J)[:on_vertices]
  push!(orbit_counts, length(orbits(Aut)))
end
