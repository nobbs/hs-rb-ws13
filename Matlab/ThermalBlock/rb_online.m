function U_rb = rb_online(Ak, F, Z, mu)
	num_Ak = size(Ak, 2);
	P = length(mu);

	A_rb = Z' * Ak{num_Ak} * Z;
	for i = 1:P
		A_rb = A_rb + mu(i) * Z' * Ak{i} * Z;
	end

	F_rb = Z' * F;

	U_rb = A_rb \ F_rb
end
