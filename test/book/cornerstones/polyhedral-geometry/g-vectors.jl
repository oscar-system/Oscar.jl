n_vertices = 30;
n_samples = 100;
g_vectors = Array{Int32}(undef,n_samples,2);

for i=1:n_samples
    RS = rand_spherical_polytope(6,n_vertices)
    g = g_vector(RS)
    g_vectors[i,1] = g[3] # notice index shift as Julia counts from 1
    g_vectors[i,2] = g[4]
end

using Plots
scatter(g_vectors[:,1], g_vectors[:,2],
        xlabel="g_2", ylabel="g_3", legend=false);

# output
Plot{Plots.GRBackend() n=1}
