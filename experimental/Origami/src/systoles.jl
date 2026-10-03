function systolic_ratio(ori::Origami)
  return Float64(GAP.Globals.SystolicRatio(GapObj(ori)).systolic_ratio)
end

function systolic_ratio_bigger_one_over_pi_in_h11(deg::Int)
  oris = Origami[]
  for d in 1:deg
    oris_deg = origamis([1, 1], d)
    append!(oris, oris_deg)
  end
  systolic_ratios = map(systolic_ratio, oris)
  threshold = 1 / pi
  amount_bigger = count(x -> x > threshold, systolic_ratios)
  total = length(oris)
  percentage = (amount_bigger / total) * 100
  rounded_percentage = round(percentage; digits=2)
  return (total, amount_bigger, "$rounded_percentage%")
end
