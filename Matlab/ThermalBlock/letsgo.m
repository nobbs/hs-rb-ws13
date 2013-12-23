% Lösen
[A, Ak, f] = assemble(p, e, t, bx, by, Mu);
U = A \ f;

% Plotten
pdesurf(p, t, U);
xlabel('x');
ylabel('y');

% Funktional "auswerten" (Mittelwert der u_ij auf dem unteren Rand... ?)
bottom_p_ind = find(p(2, :) == 0);
s = sum(U(bottom_p_ind)) / size(bottom_p_ind, 2)

