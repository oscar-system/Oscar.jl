using Plots
l(a, b) = (log(abs(a + sqrt(3)*b)), log(abs(a - sqrt(3)*b)))
pts = filter!(z -> maximum(abs.(z)) < 3,
              [l(a, b) for a in -50:50 for b in -50:50])
scatter(pts, framestyle=:origin, leg = false, mc = :blue, alpha = 0.25)
x = range(-3, 3, 100); y = -x; plot!(x, y, color = :black)
units = [ i .* l(2, 1) for i in -2:2]
scatter!(units, color = :red)
# output
Plot{Plots.GRBackend() n=3}
