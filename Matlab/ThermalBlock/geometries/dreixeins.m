global tgl_plot;
global tgl_test;
global tgl_pause;

%% Domain laden
load('geom_dreixeins.mat');
bx = 3;
by = 1;

n_train = 1e4;
n_test = 1e3;

tgl_plot = 0;
tgl_test = 0;
tgl_pause = 0;
