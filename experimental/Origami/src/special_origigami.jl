function random_origami(d::Integer)
    symm_group = symmetric_group(d)
    h = rand(symm_group)
    v = rand(symm_group)
    perm_group = permutation_group(d, [h, v])
    while !(transitivity(perm_group) > 0)
        v = rand(symm_group)
    end
    return origami_disconnected(h, v, d)
end

function staircase_origami(length::Integer, height::Integer, steps::Integer)
    h = @perm ()
    v = @perm ()
    sum_l_h = length + height
    G = symmetric_group(steps * sum_l_h)

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
