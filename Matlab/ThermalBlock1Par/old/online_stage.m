Mu_test = generate_parameter_grid(mu_min, mu_max, 1, 1, 'mc_lin');
[Ak, F] = assemble(p, e, t, bx, by);

A = Ak(:, :, end);
X = sum(Ak, 3);
for i = 1:1
	A = A + Mu_test(i) * Ak(:, :, i);
end

U_rb = rb_online(Mu_test, Z, Ak, F);
U_rb_d = assemble_sol(U_rb);
S_rb = F' * U_rb

U_fe = A \ F;
U_fe_d = assemble_sol(U_fe);
S_fe = F' * U_fe

pdesurf(p, t, abs(U_rb_d - U_fe_d));
xlabel('x');
ylabel('y');
colormap jet;
