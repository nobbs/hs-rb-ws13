% Grid = 0:0.01:1
% ep = 0.1

% FE-Approx
[Ak, F] =assemble(Grid);
A = ep * Ak(:,:,1) + Ak(:,:,2);
U_fe = A \ F;
figure();

hold on;
UU_fe = [0; U_fe; 0];
plot(Grid, U_fe, 'r');

% exakte Lösung
U_ex = exakteLoesung(Grid, ep);
plot(Grid, U_ex, 'g');

hold off;
figure();
plot(Grid, U_fe - U_ex);
