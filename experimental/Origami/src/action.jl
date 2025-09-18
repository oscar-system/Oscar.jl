# see https://arxiv.org/pdf/1809.10327.pdf Lemma 5.4

function action_s(o::Origami)
    h = horizontal_perm(o)
    v = vertical_perm(o)
    return origami_disconnected(v^-1, h, degree(o))
end

function action_t(o::Origami)
    h = horizontal_perm(o)
    return origami_disconnected(h, h^-1 * vertical_perm(o), degree(o))
end

function action_t_inv(o::Origami)
    h = horizontal_perm(o)
    v = vertical_perm(o)
    return origami_disconnected(h, h * v, degree(o))
end

function action_s_inv(o::Origami)
    h = horizontal_perm(o)
    v = vertical_perm(o)
    return origami_disconnected(v, h^-1, degree(o))
end

function action_sl2(A::ZZMatrix, o::Origami)
    # TODO necessary? loses performance
    sl2z = special_linear_group(2, ZZ)
    if !(A in sl2z)
        throw("Matrix has to be an element of SL2(Z)!")
    end

    # TODO implement this in Oscar when implementing the modular subgroup pkg
    # because this is really slow and ugly
    gap_origami = GAP.Globals.ActionOfSL2(GapObj(A), GapObj(o))
    d = GAP.Globals.DegreeOrigami(gap_origami)::Int
    G = symmetric_group(d)
    h = PermGroupElem(G, GAP.Globals.HorizontalPerm(gap_origami))
    v = PermGroupElem(G, GAP.Globals.VerticalPerm(gap_origami))
    return origami_disconnected(h, v, degree(o))
end
