using Pkg
Pkg.add("Tally"; io=devnull)


# pre-run some computation so that the same code doesn't affect the seeding when
# running the proper booktests
Qx, x = QQ["x"];
K, a = number_field(x^2 - 235, "a")
OK = ring_of_integers(K);
class_group(OK);
