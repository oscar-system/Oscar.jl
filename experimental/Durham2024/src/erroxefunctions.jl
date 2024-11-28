function embedding(C::EffectiveCartierDivisor)
	I=ideal_sheaf(C);
	C1,inc_C = sub(I);
	return inc_C
end

function reduce_to_surface()
end