function intersect(TV1::TropicalVarietySupertype{M, EMB}, TV2::TropicalVarietySupertype{M, EMB}) where {M, EMB}
    pm_tv1 = pm_object(TV1)
    pm_tv2 = pm_object(TV2)
    pm_cr = Polymake.fan.PolyhedralComplex(Polymake.fan.common_refinement(pm_tv1, pm_tv2))
    return TropicalVariety{M, EMB}(pm_cr)
end
