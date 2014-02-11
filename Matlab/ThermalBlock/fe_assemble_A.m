%% fe_assemble_A:
function A = fe_assemble_A(Ak, mu)
	% Mit 1 auffüllen, falls zu kurz...
	num_Ak = length(Ak);
	num_mu = length(mu);
	if num_mu < num_Ak
		mu_tmp = ones(num_Ak, 1);
		mu_tmp(1:num_mu) = mu;
		mu = mu_tmp;
	end
	% aufaddieren
	A = mu(1) * Ak{1};
	for j = 2:num_Ak
		A = A + mu(j) * Ak{j};
	end
end
