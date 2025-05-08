# "G-sets of permutation groups"


  # natural constructions (determined by the types of the seeds)
  G = symmetric_group(10)

  Omega = gset(G)
  @benchmark order(stabilizer(Omega, 1)[1])
  @benchmark order(stabilizer(Omega, Set([1, 2]))[1])
  @benchmark order(stabilizer(Omega, [1, 2])[1])
  @benchmark order(stabilizer(Omega, (1, 2))[1])

  Omega = gset(G, [Set([1, 2])])  # action on unordered pairs
  @benchmark order(stabilizer(Omega, Set([1, 2]))[1])
  @benchmark order(stabilizer(Omega, Set([Set([1, 2]), Set([1, 3])]))[1])
  @benchmark order(stabilizer(Omega, [Set([1, 2]), Set([1, 3])])[1])
  @benchmark order(stabilizer(Omega, (Set([1, 2]), Set([1, 3])))[1])


  Omega = gset(G, [[1, 2]])  # action on ordered pairs
  @benchmark order(stabilizer(Omega, [1, 2])[1])
  @benchmark order(stabilizer(Omega, Set([[1, 2], [1, 3]]))[1])
  @benchmark order(stabilizer(Omega, [[1, 2], [1, 3]])[1])
  @benchmark order(stabilizer(Omega, ([1, 2], [1, 3]))[1])


  Omega = gset(G, [(1, 2)])  # action on ordered pairs (repres. by tuples)
  @benchmark order(stabilizer(Omega, (1, 2))[1])
  @benchmark order(stabilizer(Omega, Set([(1, 2), (1, 3)]))[1])
  @benchmark order(stabilizer(Omega, [(1, 2), (1, 3)])[1])
  @benchmark order(stabilizer(Omega, ((1, 2), (1, 3)))[1])


  # constructions by explicit action functions
  G = symmetric_group(6)
  omega = [0,1,0,1,0,1]
  Omega = gset(G, permuted, [omega, [1,2,3,4,5,6]])

  @benchmark order(stabilizer(Omega, omega)[1])
  @benchmark order(stabilizer(Omega, Set([omega, [1,0,0,1,0,1]]))[1])
  @benchmark order(stabilizer(Omega, [omega, [1,0,0,1,0,1]])[1])
  @benchmark order(stabilizer(Omega, (omega, [1,0,0,1,0,1]))[1])
