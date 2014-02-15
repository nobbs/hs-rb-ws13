% Path erweitern
addpath('./geometries');
addpath('./tools');

% Gespeicherten Workspace laden
clear;
load('3par_workspace/matlab.mat');
print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N);

% plotten
semilogy(tab_N, tab_Delta_s_N_max, 's-', 'MarkerSize', 4);
title('Entwicklung von \Delta^s_{N,max} in Abhängigkeit von N');
ylabel('\Delta^s_{N,max}');
xlabel('N');
