#@testset "julie tests for partitions" begin
#  @test partitions(3) isa Vector
#end

@testset "combinatorics/partitions.jl" begin

	# Check some stupid cases
	@test num_partitions(0) == 1
	@test num_partitions(1) == 1
	@test num_partitions(0,0) == 1
	@test num_partitions(1,0) == 0
	@test num_partitions(1,1) == 1
	@test num_partitions(0,1) == 0


	@test num_partitions(ZZ(991)) == ZZ(16839773100833956878604913215477)

	@test num_partitions(ZZ(1991),ZZ(170)) == ZZ(22381599503916828837298114953756766080813312)
	@test num_partitions(ZZ(1991),ZZ(1000)) == ZZ(16839773100833956878604913215477)
	@test num_partitions(ZZ(1991),ZZ(670)) == ZZ(3329965216307826492368402165868892548)
	@test num_partitions(ZZ(1991),ZZ(1991)) == ZZ(1)
	@test num_partitions(ZZ(1991),ZZ(1)) == ZZ(1)

	# Constructors
	@test Partition(2,2,1) == Partition([2,2,1])
	@test Partition(1) == Partition([1])
	@test Partition() == Partition([])

	@test Partition{Int8}(2,2,1) == Partition(Int8[2,2,1])
	@test typeof(Partition{Int8}(2,2,1)) == typeof(Partition(Int8[2,2,1]))
	@test Partition{Int8}(1) == Partition(Int8[1])
	@test typeof(Partition{Int8}(1)) == typeof(Partition(Int8[1]))
	@test Partition{Int8}() == Partition(Int8[])
	@test typeof(Partition{Int8}()) == typeof(Partition(Int8[]))

	# Unrestricted partitions
	check = true
	N = 0:20
	for n in N
		P = partitions(n)
		# check that number of partitions is correct
		if length(P) != num_partitions(n)
			check = false
			break
		end
		# check that all partitions are distinct
		if P != unique(P)
			check = false
			break
		end
		# check that partitions are really partitions of n
		for lambda in P
			if sum(lambda) != n
				check = false
				break
			end
		end
	end
	@test check==true

	# ascending partitions, similar test as above
	check = true
	N = 0:20
	for n in N
		# go through all algorithms
		for a in [ "ks", "m" ]
			P = ascending_partitions(n,alg=a)
			# check that number of partitions is correct
			if length(P) != num_partitions(n)
				check = false
				break
			end
			# check that all partitions are distinct
			if P != unique(P)
				check = false
				break
			end
			# check that partitions are really partitions of n
			for lambda in P
				if sum(lambda) != n
					check = false
					break
				end
			end
		end
		if check==false
			break
		end
	end
	@test check==true

	# k-restricted partitions
	check = true
	N = 0:20
	for n in N
		for k = 0:n+1
			P = partitions(n,k)
			# check that partitions have k parts
			if length(P) !=0 && unique([ length(p) for p in P ]) != [k]
				check = false
				break
			end
			# check that number of partitions is correct
			if length(P) != num_partitions(n,k)
				check = false
				break
			end
			# check that all partitions are distinct
			if P != unique(P)
				check = false
				break
			end
			# check that partitions are really partitions of n
			for lambda in P
				if sum(lambda) != n
					check = false
					break
				end
			end
		end
		if check==false
			break
		end
	end
	@test check==true

	# k-restricted partitions with lower and upper bounds
	check = true
	N = 0:20
	for n in N
		for k = 0:n+1
			for l1 = 0:n
				for l2 = l1:n
					P = partitions(n,k,l1,l2)
					# check that partitions have k parts
					if length(P) !=0 && unique([ length(p) for p in P ]) != [k]
						check = false
						break
					end
					# check that all partitions are distinct
					if P != unique(P)
						check = false
						break
					end
					# check that partititons are really partitions of n
					for lambda in P
						if sum(lambda) != n
							check = false
							break
						end
					end
					#check that all parts are inside the bounds
					for lambda in P
						for part in lambda
							if part>l2 || part<l1
								check = false
								break
							end
						end
					end
				end
			end
		end
		if check==false
			break
		end
	end
	@test check==true

	# k-restricted partitions with lower and upper bounds and distinct parts
	check = true
	N = 0:20
	for n in N
		for k = 0:n+1
			for l1 = 0:n
				for l2 = l1:n
					P = partitions(n,k,l1,l2,z=1)
					# check that partitions have k parts
					if length(P) !=0 && unique([ length(p) for p in P ]) != [k]
						check = false
						break
					end
					# check that all partitions are distinct
					if P != unique(P)
						check = false
						break
					end
					# check that partititons are really partitions of n
					for lambda in P
						if sum(lambda) != n
							check = false
							break
						end
					end
					#check that all parts are inside the bounds
					for lambda in P
						for part in lambda
							if part>l2 || part<l1
								check = false
								break
							end
						end
					end
					#check that all parts are distinct
					for lambda in P
						if lambda != unique(lambda)
							check = false
							break
						end
					end
				end
			end
		end
		if check==false
			break
		end
	end
	@test check==true

	#v-mu-restricted partitions
	check = true
	N = 0:20
	for n in N
		for k = 1:n+1
			for i = 1:20
				v = convert(Array{Integer,1}, rand(1:20,i))
				unique!(v)
				sort!(v)
				mu = convert(Array{Integer,1}, rand((0:k).+1, length(v)))
				P = partitions(mu, n, v, k)
				# check that partitions have k parts
				if length(P) !=0 && unique([ length(p) for p in P ]) != [k]
					check = false
					break
				end
				# check that all partitions are distinct
				if P != unique(P)
					check = false
					break
				end
				# check that partitions are really partitions of n
				for λ in P
					if sum(λ) != n
						check = false
						break
					end
				end
				#=
				# check that partitions are constrained by v and mu
				for λ in P
					partcount = partition_to_partcount(λ)
					for j = 1:length(partcount)
						if partcount[j] != 0
							ind = findfirst(isequal(j),v)
							if ind==nothing || mu[ind] < partcount[j]
								check = false
								break
							end
						end
					end
				end =#
			end
		end
		if check==false
			break
		end
	end
	@test check==true
	#check for a few examples if the correct v-mu-restricted Partitions are found
	check = true
	N = 1:20
	K = 1:20
	for n in N
		for k in K
			if partitions(n,k) != partitions(fill(n,n),n,collect(1:n),k)
				check = false
			end
			if partitions(n, k, 1, n, z=1) != partitions(fill(1,n),n,collect(1:n),k)
				check = false
			end
			for l1 = 1:n
				for l2 = l1:n
					if partitions(n, k, l1, l2, z=1) != partitions(fill(1,l2-l1+1),n,collect(l1:l2),k)
						check = false
					end
					if partitions(n, k, l1, l2) != partitions(fill(n,l2-l1+1),n,collect(l1:l2),k)
						check = false
					end
				end
			end
		end
	end
	@test_throws ArgumentError partitions([1,2], 5, [1,2,3], 2)
	@test_throws ArgumentError partitions([1,2,3], -1, [1,2,3], 2)
	@test_throws ArgumentError partitions([1,2,3], 5, [1,2], 0)

	# Dominance order
	@test dominates(Partition([4,2]), Partition([3,2,1])) == true
	@test dominates(Partition([4,1,1]), Partition([3,3])) == false
	@test dominates(Partition([3,3]), Partition([4,1,1])) == false
	@test dominates(Partition([5]), Partition([2,2,2])) == false

	# Conjugate partition
	@test conjugate(Partition([6,4,3,1])) == Partition([4, 3, 3, 2, 1, 1])
	@test conjugate(Partition([5,4,1])) == Partition([3, 2, 2, 2, 1])
	@test conjugate(Partition([])) == Partition([])
	check = true
	for p in partitions(10)
		if conjugate(conjugate(p)) != p
			check = false
			break
		end
	end
	@test check==true

	#=
	# partcount_to_partition
	@test partcount_to_partition([2,0,1]) == Partition([3,1,1])
	@test partcount_to_partition(Int[]) == Partition([])
	@test partcount_to_partition([0]) == Partition([])
	@test partcount_to_partition([3,2,1]) == Partition([3,2,2,1,1,1])

	# partition_to_partcount
	@test partition_to_partcount(Partition([5,3,3,3,2,1,1])) == [2,1,3,0,1]
	@test partition_to_partcount(Partition([])) == []
	@test partition_to_partcount(Partition([3,2,1])) == [1,1,1]

	
	check = true
	N = 0:20
	for n in N
		for λ in partitions(n)
			if λ != partcount_to_partition(partition_to_partcount(λ))
				check=false
				break
			end
		end
		if check==false
			break
		end
	end
	@test check==true =#

end

