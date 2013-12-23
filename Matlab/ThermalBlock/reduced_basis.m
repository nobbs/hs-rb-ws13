function [Z] = reduced_basis(p, e, t, bx, by, Xi_train, n_p, Mu_bar)
tic;
n = n_p;
N_max = 100;
tol = 1e-6;

Snapshots = zeros(n, 1);
Z = zeros(n, 1);
Choosen = []

% Ersten Snapshot berechnen
disp(['snapshot with mu'])
Xi_train(:, 1)
[A, Ak, F, X] = assemble(p, e, t, bx, by, Xi_train(:, 1));
Snapshots(:, 1) = A \ F;
Z(:, 1) = Snapshots(:, 1) / (sqrt(Snapshots(:, 1).' * X * Snapshots(:, 1)));
Choosen(1) = 1;
MaxErrors = [];
toc

for i = 2:N_max
	tic;
	N = i - 1;
	P = 8;
	id = 1;
	err = 0;

	H = zeros(n, 1 + (P+1) * N);
	H(:, 1) = F;
	for j = 1:P+1
		H(:, 1+j:P+1:end) = - Ak(:, :, j) * Z;
	end
	G = H.' * X * H;

	Errors = zeros(size(Xi_train, 2), 1);
	for j = 1 : size(Xi_train, 2)
		if any(Choosen == j)
			continue;
		end

		mu_check = Xi_train(:, j);

		% Lösen
		A_rb_check = zeros(size(Z, 2));
		F_rb_check = Z.' * F;
		for k = 1:P
			A_rb_check = A_rb_check + mu_check(k) * Z.' * Ak(:, :, k) * Z;
		end
		A_rb_check = A_rb_check + Z.' * Ak(:, :, P + 1) * Z;
		U_rb_check = A_rb_check \ F_rb_check;
		% pdesurf(p, t, Z * U_rb_check)
		% waitforbuttonpress;

		% Fehler berechnen
		Eps = zeros(1 + (P+1) * N, 1);
		Eps(1) = 1;
		for k = 1:P
			Eps(1+k:P+1:end) = mu_check(k) * U_rb_check;
		end
		Eps(1+P:P+1:end) = U_rb_check;

		% Fehler
		resid = sqrt(Eps.' * G * Eps);
		alpha_LB = min([mu_check ./ Mu_bar]);
		Err_en = resid / sqrt(alpha_LB);
		Errors(j) = Err_en;
		if Err_en > err
			err = Err_en;
			id = j;
		end
	end

	if err < tol
		break;
	end

	MaxErrors(end + 1) = err;
	Choosen(end + 1) = id;
	% plot(2:i, MaxErrors);

	disp(['New Snapshot for N = ', num2str(N + 1), ', max_error = ', num2str(err) ,' for mu = '])
	Xi_train(:, id)
	[A, Ak, F, X] = assemble(p, e, t, bx, by, Xi_train(:, id));
	Snapshots(:, i) = A \ F;

	z = Snapshots(:, i);
	for j = 1:(i-1)
		z = z - Snapshots(:, i).' * X * Z(:, j) * Z(:, j);
	end
	Z(:, i) = z / (sqrt(z.' * X * z));
	toc
end


% N_max = size(Xi_train, 2);
% Snapshots = zeros(n_p, N_max);
% Snapshots berechnen
% for i = 1:N_max
	% disp('snapshot')
	% [A, Ak, F, X] = assemble(p, e, t, bx, by, Xi_train(:, i));
	% Snapshots(:, i) = A \ F;
% end

% Gram-Schmidt-Orthonormalisierung
% Z = zeros(n_p, N_max);
% Z(:, 1) = Snapshots(:, 1) / (sqrt(Snapshots(:, 1)' * X * Snapshots(:, 1)));
% for i = 2:N_max
	% z = Snapshots(:, i);
	% for j = 1:(i-1)
		% z = z - (Snapshots(:, i)' * X * Z(:, j)) * Z(:, j);
	% end
	% Z(:, i) = z / (sqrt(z' * X * z));
% end

end
