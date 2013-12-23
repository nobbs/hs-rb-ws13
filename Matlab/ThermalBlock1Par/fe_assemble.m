% Assembliert die Steifigkeitsmatrix.
%
% - Neumann-Bedingungen verändern in diesem Fall nur die rechte Seite F.
% - Dirichlet-Bedingung wird berücksichtigt, indem die betroffenen Punkte aus der
%   Steifigkeitsmatrix gestrichen werden (dadurch wird diese wieder symmetrisch)
%   und kann nach dem Berechnen der Lösung schnell wieder hinzugefügt werden
%   (via U = B * Uo + Ud, wobei Uo die Lösung ohne Randpunkte ist)
function [Ak, F, B, Ud] =  fe_assemble(P, E, T, bx, by)
    num_points = size(P, 2);        % Anzahl Gitterpunkte
    num_blocks = bx * by;           % Anzahl Blöcke
    num_points_dbc = 0;             % Anzahl Gitterpunkte, von Dirichlet-RB betroffen

    Ak = cell(num_blocks, 1);
    for j = 1:num_blocks
        Ak{j} = sparse(num_points, num_points);
    end
    F = sparse(num_points, 1);

    % Steifigkeitsmatrizen aufbauen
    for i = 1:size(T, 2)
        % Subdomain legt fest, in welcher Matrix der Eintrag landet
        subdomain = T(4, i);
        % Elementsteifigkeitsmatrix aufbauen
        [Sk, fk] = fe_assemble_element([P(:, T(1, i)), P(:, T(2, i)), P(:, T(3, i))]);
        % Elementsteifigkeitsmatrix einbauen
        Ak{subdomain}(T(1:3, i), T(1:3, i)) = Ak{subdomain}(T(1:3, i), T(1:3, i)) + Sk;
        F(T(1:3, i)) = F(T(1:3, i)) + fk;
    end

    % Von Dirichlet-RB betroffene Punkte suchen
    P_dirichlet = [];
    for i = find(P(2, :) == 1)
        num_points_dbc = num_points_dbc + 1;
        P_dirichlet(end + 1) = i;
    end

    num_points_no_dbc = num_points - num_points_dbc;
    range = 1:num_points;
    range(P_dirichlet) = [];

    % Betroffene Punkte streichen
    for j = 1:num_blocks
        Ak{j} = Ak{j}(range, range);
    end
    F = F(range);

    Ud = sparse(num_points, 1);
    B = sparse(num_points, num_points_no_dbc);
    B(range, :) = eye(num_points_no_dbc, num_points_no_dbc);
end

function [Sk, fk] = fe_assemble_element(C)
% Liegt der Punkt auf der unteren Kante? -> F bekommt Einträge
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
