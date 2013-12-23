function [G] = generate_parameter_grid(mu_min, mu_max, P, m, method)
	if method == 'mc_lin'
		G = mu_min * ones(P, m) + (mu_max - mu_min) * rand(P, m) .* ones(P, m);
	elseif method == 'mc_log'
		mu_log_min = log(mu_min);
		mu_log_max = log(mu_max);
		G = mu_log_min * ones(P, m) + (mu_log_max - mu_log_min) * rand(P, m) .* ones(P, m);
	end
end
