function [Ak, F] = assemble(Grid)
n = length(Grid)
Ak = zeros(n, n, 2);
F = zeros(n, 1);

% A^1, A^2 und F basteln
for i = 2:(n-1)
	Ak(i, i, 1) = 1 / (Grid(i) - Grid(i - 1)) + 1 / (Grid(i + 1) - Grid(i));
	Ak(i, i + 1, 1) = - 1 / (Grid(i + 1) - Grid(i));
	Ak(i, i - 1, 1) = - 1 / (Grid(i) - Grid(i - 1));
	F(i) = 0.5 * (Grid(i + 1) - Grid(i - 1));
end
Ak(:, :, 2) = diag(-0.5 * ones(n - 1, 1), -1) + diag(0.5 * ones(n - 1, 1), 1);
Ak(1, :, :) = 0;
Ak(n, :, :) = 0;
Ak(1, 1, :) = 1;
Ak(n, n, :) = 1;

end
