global tgl_plot;
global tgl_test;
global tgl_pause;
global tgl_train_eq_test;
global tgl_logplot;

%% Domain laden
load('geom_zweixzwei.mat');
bx = 2;
by = 2;

n_train = 1e5;
n_test = 1e5;

tgl_plot = 0;
tgl_test = 1;
tgl_pause = 0;
tgl_train_eq_test = 1
tgl_logplot = 0
