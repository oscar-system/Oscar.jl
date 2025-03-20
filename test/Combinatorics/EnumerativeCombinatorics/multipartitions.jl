@testset "Multipartitions" begin
	#constructors
	@test Multipartition([[3,2],[1]]) == Multipartition([Partition([3,2]),Partition([1])])
	@test Multipartition([[3,2],[]]) == Multipartition([Partition([3,2]),Partition([])])

	# multi-partitions
	check = true
	N = 0:10
	R = 1:5
	for n in N
		for r in R
			MP = multipartitions(n,r)
			# check that all multipartitions are distinct
			if MP != unique(MP)
				check = false
				break
			end
			# check that multipartititons are really multipartitions of n
			for mp in MP
				if sum(mp) != n
					check = false
					break
				end
			end
			# check that all multisetpartitions have k parts
			if length(MP) !=0 && unique([ length(mp) for mp in MP ]) != [r]
				check = false
				break
			end
		end
	end
	@test check==true

	#counting
	for n=0:5
		for k=1:n+1
			@test length(multipartitions(n,k)) == num_multipartitions(n,k)
		end
	end

end
