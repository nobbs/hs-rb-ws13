%% print_test_res: function description
function print_test_res(tab_N, tab_Delta_s_N_max, tab_eta_s_N_ave, tab_eta_s_N_max, tab_rho_S_err_N)
    disp(' ');
    fprintf('%s\t%s\t%s\t%s\t%s\n', 'N', 'Delta_s_N_max', 'eta_s_N_ave', 'eta_s_N_max', 'rho_S_err_N');
    for i = 1:length(tab_N)
        fprintf('%d\t%10.6e\t%10.6f\t%10.6f\t%10.6f\n', tab_N(i), tab_Delta_s_N_max(i), tab_eta_s_N_ave(i), tab_eta_s_N_max(i), tab_rho_S_err_N(i));
    end
    disp(' ');
end
