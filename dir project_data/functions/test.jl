n = 5
M = []
graphs = []
partitions = lvl2_leaf_partition(n)
for part in partitions
    push!(M, [hybrid_vertex_mod_sym(part), part])
end

counter = 0
types = [1]
for m in M
    for (i, rec_vert) in enumerate(m[1])
        print(rec_vert, m[2], "\n")
        G = draw_network(rec_vert, m[2])
        name = string.(n)*"_"*string(m[2])*"_ideal_"*string(i)
        push!(graphs, [G, name])
        counter = counter + 1
    end
    push!(types, counter)
end

Oscar.save("dir_project_data/network_data/"*string(n)*"_leaves_types_data", types)

results = compare_networks(graphs)

Oscar.save("dir_project_data/network_data/"*string(n)*"_leaves_data", results)

display(results)
print(types[4])