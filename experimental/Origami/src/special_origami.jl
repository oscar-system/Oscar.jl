function random_origami(d::Integer)
    symm_group = symmetric_group(d)
    h = rand(symm_group)
    v = rand(symm_group)
    perm_group = permutation_group(d, [h, v])
    while !(transitivity(perm_group) > 0)
        v = rand(symm_group)
        perm_group = permutation_group(d, [h, v])
    end
    return origami_disconnected(h, v, d)
end

function staircase_origami(length::Integer, height::Integer, steps::Integer)
    sum_l_h = length + height
    G = symmetric_group(steps * sum_l_h)
    h = one(G)
    v = one(G)

    for step in 1:(steps - 1)
        h_start = (step - 1) * (sum_l_h) + 1
        h_end = (step - 1) * (sum_l_h) + length
        h = h * cperm(G, h_start:h_end)

        v_start = step * length + (step - 1) * height
        v_end = step * (sum_l_h) + 1
        v = v * cperm(G, v_start:v_end)
    end

    h_start = (steps - 1) * (sum_l_h) + 1
    h_end = (steps - 1) * (sum_l_h) + length
    h = h * cperm(G, h_start:h_end)

    v_start = steps * length + (steps - 1) * height
    v_end = steps * sum_l_h
    v = v * cperm(G, v_start:v_end)

    return normal_form(origami(h, v))
end

function x_origami(half_degree::Integer)
    G = symmetric_group(2*half_degree)
    sigma_h = cperm(G, 1:2*half_degree)
    sigma_v = one(G)

    for tile in 1:half_degree
        sigma_v = sigma_v * cperm(G, [2 * tile - 1, 2 * tile])
    end

    return normal_form(origami(sigma_h, sigma_v))
end

function elevator_origami(length::Integer, height::Integer, steps::Integer)
    sigma_h = cperm()
    sigma_h_step = Vector{PermGroupElem}(undef, steps)

    for step in 1:steps
        sigma_h_step[step] = cperm((step-1)*(length+height)+1:(step-1)*(length+height)+length)
        sigma_h = sigma_h * sigma_h_step[step]
    end

    sigma_v = cperm()
    sigma_v_step = Vector{PermGroupElem}(undef, steps-1)
    for step in 1:(steps-1)
        sigma_v_step[step] = cperm(step*length+(step-1)*height:step*(length+height)+1)
        sigma_v = sigma_v * sigma_v_step[step]
    end
    last_connection = collect(steps*length+(steps-1)*height:steps*(length+height))
    push!(last_connection, 1)
    sigma_v = sigma_v * cperm(last_connection)

    return normal_form(origami(sigma_h, sigma_v))
end
