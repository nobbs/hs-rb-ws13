%% assemble_fe_A: Kombiniert die parameterunabhängigen FE-Matrizen zu der FE-
%% Steifigkeitsmatrix zum Parameter mu (in Form von Theta).
function A = assemble_fe_A(Ak, Theta)
    Q_a = size(Ak, 1);
	A = Theta(1) * Ak{1};
	for j = 2:Q_a
		A = A + Theta(j) * Ak{j};
	end
end
