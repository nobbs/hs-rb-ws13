global tgl_plot;
global tgl_test;
global tgl_pause;
global tgl_train_eq_test;
global tgl_logplot;
global tgl_print;

%% Domain laden
load('geom_zweixeins.mat');
bx = 2;
by = 1;

n_train = 1e4;
n_test = 1e4;

tgl_plot = 1
tgl_test = 1
tgl_pause = 1
tgl_train_eq_test = 1
tgl_logplot = 1
tgl_print = 0
