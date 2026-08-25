function draw_path(p)
    return [[p[i], p[i+1]] for i in 1:length(p)-1]
end

function add_leaves(sum, path)
    path = path[2: length(path) - 1]
    edges = []
    for i in 1:length(path)
        push!(edges, [path[i], sum + i])
    end
    return edges
end

function draw_network(ret, par)
    edges = Vector{Vector{Int}}()
    v_1 = ret[1]
    v_2 = ret[2]

    if v_1 > v_2
        temp = v_2
        v_2 = v_1
        v_1 = temp
    end

    p_1 = [i for i in 1:par[1] + 2]
    p_2 = vcat([1], [i for i in par[1] + 3: par[1] + par[2] + 2], [par[1] + 2]) 
    p_3 = vcat([1], [i for i in par[1] + par[2] + 3: par[1] + par[2] + par[3]+ 2], [par[1] + 2])
    paths = [p_1, p_2, p_3]

    if v_1 == 1
        q_1, q_2, q_3 = [], [], []
        if v_2 in p_1
            q_1 = p_3
            q_2 = p_2
            q_3 = p_1
        end
        if v_2 in p_2
            q_1 = p_1
            q_2 = p_3
            q_3 = p_2
        end
        if v_2 in p_3
            q_1 = p_1
            q_2 = p_2
            q_3 = p_3
        end
        append!(edges, draw_path(q_1))
        append!(edges, draw_path(q_2))
        i = findfirst(==(v_2), q_3)
        A1 = vcat(q_3[1:i-1], [v_2])
        A2 = q_3[i:end]
        append!(edges, draw_path(A1))
        append!(edges, draw_path(reverse(A2)))
    else
        if (v_1 in p_1) && (v_2 in p_2)
            q_1 = p_3
            q_2 = p_1
            q_3 = p_2
        end

        if (v_1 in p_1) && (v_2 in p_3)
            q_1 = p_2
            q_2 = p_1
            q_3 = p_3
        end

        if (v_1 in p_2) && (v_2 in p_3)
            q_1 = p_1
            q_2 = p_2
            q_3 = p_3
        end

        append!(edges, draw_path(q_1))

        i = findfirst(==(v_1), q_2)
        A1 = vcat(q_2[1:i-1], [v_1])
        A2 = q_2[i:end]
        append!(edges, draw_path(A1))
        append!(edges, draw_path(reverse(A2)))

        j = findfirst(==(v_2), q_3)
        B1 = vcat(q_3[1:j-1], [v_2])
        B2 = q_3[j:end]
        append!(edges, draw_path(B1))
        append!(edges, draw_path(reverse(B2)))
    end

    sum = par[1] + par[2] + par[3] + 2

    for i in 1:3
        append!(edges, add_leaves(sum, paths[i]))
        sum = sum + par[i]
    end
    g = graph_from_edges(Directed, edges)
    n = phylogenetic_network(g)
    M = cavender_farris_neyman_model(n)

    return M
end

draw_network([5,6], [2,1,1])