% Grid (im kritischen Bereich genauer aufgelöst)
% Grid = unique([0:0.1:1]);
Grid = unique([0:0.01:0.9, 0.9:0.001:1]);
n_with_borders = length(Grid);
n = n_with_borders;

% Parameter-Daten
mu_min = 0.01;
mu_max = 0.5;
mu_bar = 1;

D_box = [mu_min, mu_max];
D_log_box = [log(mu_min), log(mu_max)];

% äquidistantes Parameter-Gitter:
steps = 1000;
% G_lin_aq = mu_min + (0:steps) * (mu_max - mu_min) / steps;
% G_log_mc = mu_min * exp(log(mu_max / mu_min) * (0:steps) / steps);
% G_lin_mc = sort(unique(mu_min * ones(1, steps) + (mu_max - mu_min) * rand(1, steps)));
G_log_mc = sort(unique(mu_min * exp(log(mu_max / mu_min) * rand(1, steps))));

% Offline-Stage ausführen
Xi_train = G_log_mc;

[Ak, F] = assemble(Grid);
X = mu_bar * Ak(:, :, 1) + Ak(:, :, 2);

N_max = 50;
tolerance = 1e-4;
Choosen = [];

Snapshots = zeros(n, 1);
Z = zeros(n, 1);

% Ersten Snapshot berechnen
disp(['snapshot with mu = ', num2str(Xi_train(floor(steps / 2)))])
A = Xi_train(floor(steps / 2)) * Ak(:, :, 1) + Ak(:, :, 2);
Choosen(1) = Xi_train(floor(steps / 2));
Snapshots(:, 1) = A \ F;
Z(:, 1) = Snapshots(:, 1) / (sqrt(Snapshots(:, 1)' * X * Snapshots(:, 1)));
Choosen(1) = floor(steps / 2);

for i = 2:N_max
	N = i - 1;
	% disp('Searching biggest error');
	id = 1;
	err_en = 0;
	err_s = 0;
	U_rb_N = 0;
	S_N = 0;

	Errors = zeros(length(Xi_train), 1);
	for j = 1 : length(Xi_train)
		if any(Choosen == j)
			continue;
		end

		H = zeros(n, 1 + 2 * N);
		H(:, 1) = F;
		H(:, 2:2:end) = - Ak(:, :, 1) * Z;
		H(:, 3:2:end) = - Ak(:, :, 2) * Z;
		G = H.' * (X^-1).' * H;

		mu_check = Xi_train(j);

		% Lösen
		A_rb_check = Z.' * (mu_check * Ak(:, :, 1) + Ak(:, :, 2)) * Z;
		F_rb_check = Z.' * F;
		U_rb_check = A_rb_check \ F_rb_check;

		% u_n_norm = sqrt(U_rb_check.' * (Z.' * X * Z) * U_rb_check)

		% Fehler berechnen
		Eps = zeros(1 + 2 * N, 1);
		Eps(1) = 1;
		Eps(2:2:end) = mu_check * U_rb_check;
		Eps(3:2:end) = U_rb_check;

		% Fehler
		resid = real(sqrt(Eps.' * G * Eps));
		alpha_LB = min([mu_check / mu_bar, 1 / mu_bar]);
		Err_en = resid / sqrt(alpha_LB);
		Err_s = Err_en^2;
		Errors(j) = Err_en;
		if Err_en > err_en
			U_rb_N = Z * U_rb_check;
			S_N = F' * U_rb_N;
			err_en = Err_en;
			err_s = Err_s;
			id = j;
		end
	end

	semilogy(Xi_train, Errors);
	hold on;
	for y = Xi_train(Choosen)
		line([y y], [1e-6 100], 'LineStyle', '--', 'Color', 'r');
	end
	hold off;
	xlim([0, 0.5])
	ylim([1e-6, 100])
	waitforbuttonpress;

	if err_en < tolerance
		break;
	end

	disp(['snapshot with mu = ', num2str(Xi_train(id))])
	A = Xi_train(id) * Ak(:, :, 1) + Ak(:, :, 2);
	Snapshots(:, i) = A \ F;

	Choosen(end + 1) = id;
	disp('Fehlerschaetzer err_en')
	err_en
	disp('Echter Fehler ||u - u_N||_X')
	sqrt((Snapshots(:, i) - U_rb_N).' * X * (Snapshots(:, i) - U_rb_N))
	disp('Fehlerschaetzer err_s')
	err_s
	disp('Echter Fehler |s - s_N|')
	S = F.' * Snapshots(:, i);
	abs(S - S_N)

	if i == 1
		Z(:, 1) = Snapshots(:, 1) / (sqrt(Snapshots(:, 1)' * X * Snapshots(:, 1)));
	else
		z = Snapshots(:, i);
		for j = 1:(i-1)
			z = z - (Snapshots(:, i)' * X * Z(:, j)) * Z(:, j);
		end
		Z(:, i) = z / (sqrt(z' * X * z));
	end
end

plot(Xi_train, Errors);
hold on;
plot(Xi_train(Choosen), 0, 'r*');
hold off;
% waitforbuttonpress;
