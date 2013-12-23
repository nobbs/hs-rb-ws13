function [A, Ak, F, X] =  assemble(P, E, T, bx, by, Mu)
n = size(P, 2);
N = bx * by;

Ak = zeros(n, n, N);
F = zeros(n, 1);

% Steifigkeitsmatrizen aufbauen
for i = 1:size(T, 2)
   p1 = P(:, T(1, i));
   p2 = P(:, T(2, i));
   p3 = P(:, T(3, i));
   sd = T(4, i);
   [Sk, fk] = assembleElement([p1, p2, p3]);
   Ak(T(1:3, i), T(1:3, i), sd) = Ak(T(1:3, i), T(1:3, i), sd) + Sk;
   F(T(1:3, i)) = F(T(1:3, i)) + fk;
end

% X = zeros(n);
% for i = 1:N
% 	X = X + Ak(:, :, i);
% end
X = sum(Ak, 3);

% Dirichlet-RB obere Kante
for i = find(P(2, :) == 1)
	Ak(i, :, :) = 0;
	[r, c] = find(T(1:3, :) == i);
	for sd = T(4, c)
		Ak(i, i, sd) = 1;
	end
end

% A aufsummieren
A = zeros(n);
A = Ak(:, :, N);
for i = 1:(N - 1)
	A = A + Mu(i) * Ak(:, :, i);
end

end

function [Sk, fk] = assembleElement(C)
% Liegt der Punkt auf der unteren Kante?
at_bottom = false;
P = [];
if C(2, 1) == 0 & C(2, 2) == 0
	at_bottom = true;
	P = [1, 2];
elseif C(2, 1) == 0 & C(2, 3) == 0
	at_bottom = true;
	P = [1, 3];
elseif C(2, 2) == 0 & C(2, 3) == 0
	at_bottom = true;
	P = [2, 3];
end

% etas
e31 = C(2, 3) - C(2, 1);
e12 = C(2, 1) - C(2, 2);
% xis
x13 = C(1, 1) - C(1, 3);
x21 = C(1, 2) - C(1, 1);

d = x21 * e31 - x13 * e12;

if (d <= 0)
    error('Flächeninhalt d ist <= 0 ?!');
end

G = [-1, -1; 1, 0; 0, 1];
Phi_inv = [e31, x13; e12, x21];

% Elementsteifigkeitsmatrix
Sk = (0.5 / d) * G * (Phi_inv * Phi_inv.') * G.';

% Bestimme die rechte Seite. Diese ist nur für den unteren Rand != 0, also
% beachte auch nur die auf diesem Rand liegenden Punkte.
fk = zeros(3, 1);
if at_bottom
	[x_a, i_a] = min(C(1, P));
	[x_e, i_e] = max(C(1, P));
	% Lambda_a =^= absteigend
	% Lambda_e =^= aufsteigend
	fk(P(i_a)) = (x_e - x_a) / 2;
	fk(P(i_e)) = (x_e - x_a) / 2;
end

end
