%% generate_parameter_grid: Erzeugt ein m-elementiges Gitter für den
%% Parameterbereich [mu_min, mu_max]^P nach der gewählten Methode.
%% Mögliche Methoden:
%% 	aq_lin: linear äquidistant	(nur für P = 1 oder 2)
%% 	aq_log: logarithmisch äquidistant  (nur für P = 1 oder 2)
%% 	mc_log: zufällig linear gleichverteilt
%% 	mc_log: zufällig logarithmisch gleichverteilt
function param_grid = generate_parameter_grid(mu_min, mu_max, P, m, method)
	if method == 'aq_lin' & P == 1
		param_grid = linspace(mu_min, mu_max, m);
	elseif method == 'aq_lin' & P == 2
		m = floor(sqrt(m));
		[X, Y] = meshgrid(linspace(mu_min, mu_max, m), linspace(mu_min, mu_max, m));
		param_grid = [X(:) Y(:)]';
	elseif method == 'aq_log' & P == 1
		param_grid = mu_min * exp(log(mu_max / mu_min) * linspace(0, 1, m));
	elseif method == 'aq_log' & P == 2
		m = floor(sqrt(m));
		[X, Y] = meshgrid(mu_min * exp(log(mu_max / mu_min) * linspace(0, 1, m)), mu_min * exp(log(mu_max / mu_min) * linspace(0, 1, m)));
		param_grid = [X(:) Y(:)]';
	elseif method == 'mc_lin'
		param_grid = mu_min * ones(P, m) + (mu_max - mu_min) * rand(P, m) .* ones(P, m);
	elseif method == 'mc_log'
		param_grid = mu_min * exp(log(mu_max / mu_min) * rand(P, m) .* ones(P, m));
	end

	if P == 1
		param_grid = sort(param_grid);
	end
end
