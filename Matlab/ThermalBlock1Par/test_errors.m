%% test_errors: function description
function [Delta_s_N_max, eta_s_N_ave, eta_s_N_max, rho_S_err_N] = test_errors(Ak, F, X, Z, Xi_test, mu_bar)
	P = length(Ak);
	N = size(Z, 2);
	num_ps = size(Ak{1}, 1);
	n_test = size(Xi_test, 2);

	% Immer wieder benötigt...
	F_rb = Z' * F;
	H = sparse(num_ps, 1 + 2 * N);
	H(:, 1) = F;
	for k = 1:P
		H(:, (1+k):P:end) = - Ak{k} * Z;
	end
	G = H' * (X \ H);

	Test_errors = zeros(n_test, 1);
	Test_effectivities = zeros(n_test, 1);
	Test_S_diff = zeros(n_test, 1);
	for j = 1:n_test
		mu = Xi_test(j);
		Theta = [mu, 1];
		A_rb = Z.' * fe_assemble_A(Ak, mu) * Z;
		U_rb = A_rb \ F_rb;

		%%% Fehler berechnen
		Eps = sparse(1 + 2 * N, 1);
		Eps(1) = 1;
		for k = 1:P
			Eps((1+k):P:end) = Theta(k) * U_rb;
		end
		residue_test = real(sqrt(Eps.' * G * Eps));

		alpha_LB_test = min([mu ./ mu_bar; 1]);
		err_en_test = residue_test / sqrt(alpha_LB_test);
		err_s_test = err_en_test^2;
		Test_errors(j) = err_s_test;

		% Effektivität
		U_fe_test = fe_assemble_A(Ak, mu) \ F;

		S_fe_test = F.' * U_fe_test;
		S_rb_test = F_rb.' * U_rb;

		Test_effectivities(j) = err_s_test / (S_fe_test - S_rb_test);
		Test_S_diff(j) = (S_fe_test - S_rb_test);
	end

	% Test_errors
	% Test_effectivities
	% Test_S_diff

	Delta_s_N_max = max(Test_errors);
	eta_s_N_ave = sum(Test_effectivities) / n_test;
	eta_s_N_max = max(Test_effectivities);
	rho_S_err_N = Delta_s_N_max / max(Test_S_diff);
end
