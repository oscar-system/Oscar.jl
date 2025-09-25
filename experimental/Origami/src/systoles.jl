function systolic_ratio(ori::Origami)
    return GAP.gap_to_julia(GAP.Globals.SystolicRatio(ori.o))[:systolic_ratio]
end

function systolic_ratio_bigger_one_over_pi_in_h11(deg::Int)
    oris = Origami[]
    for d in 1:deg
        orisDeg = origamis([1,1], d)
        append!(oris, orisDeg)
    end
    systolicRatios = map(systolic_ratio, oris)
    threshold = 1 / pi
    amountBigger = count(x -> x > threshold, systolicRatios)
    total = length(oris)
    percentage = (amountBigger / total) * 100
    rounded_percentage = round(percentage, digits=2)
    return (total, amountBigger, "$rounded_percentage%")
end