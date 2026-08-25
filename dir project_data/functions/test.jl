n = 4
M = []
graphs = []
for part in partitions
    push!(M, [hybrid_vertex_mod_sym(part), part])
end

println(M)
for m in M
    for i, rec_vert in enumerate(m[1])
        print(rec_vert, m[2], "\n")
        G = draw_network(rec_vert, m[2])
        name = string.(n)*"_"*string(m[2])*"_ideal_"*string(i)
        push!(graphs, [G, name])
    end
end

results = compare_networks(graphs)

println(results)