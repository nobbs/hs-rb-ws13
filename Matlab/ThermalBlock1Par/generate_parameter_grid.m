% Erzeugt ein m-elementiges Gitter für den Parameterbereich [mu_min, mu_max]^P
% nach der gewählten Methode.
function G = generate_parameter_grid(mu_min, mu_max, P, m, method)
	% aq_* nur für P = 1
	if method == 'aq_lin' & P == 1
		G = linspace(mu_min, mu_max, m);
	elseif method == 'aq_lin' & P == 2
		m = floor(sqrt(m));
		[X, Y] = meshgrid(linspace(mu_min, mu_max, m), linspace(mu_min, mu_max, m));
		G = [X(:) Y(:)]';
	elseif method == 'aq_log' & P == 1
		G = mu_min * exp(log(mu_max / mu_min) * linspace(0, 1, m));
	elseif method == 'aq_log' & P == 2
		m = floor(sqrt(m));
		[X, Y] = meshgrid(mu_min * exp(log(mu_max / mu_min) * linspace(0, 1, m)), mu_min * exp(log(mu_max / mu_min) * linspace(0, 1, m)));
		G = [X(:) Y(:)]';
	elseif method == 'mc_lin'
		G = mu_min * ones(P, m) + (mu_max - mu_min) * rand(P, m) .* ones(P, m);
		if P == 1
			G = sort(G);
		end
	elseif method == 'mc_log'
		G = mu_min * exp(log(mu_max / mu_min) * rand(P, m) .* ones(P, m));
		if P == 1
			G = sort(G);
		end
	end
end
